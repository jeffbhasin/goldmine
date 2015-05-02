# Functions for testing enrichment between sets of genomic ranges and generating background sets

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Test enrichment between a query set and a null set
#'
#' Provide both a query set of ranges (as a GenomicRanges, data.frame, or data.table) and a null set of ranges (same format options). Counts for each feature will be computed in both the query and null sets, and tested for significance of difference using binom.test().
#' @param query A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param null A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param features A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based. Additionally, there must be a column named "name" which will be used as a factor to divide the ranges into subsets. Each subset will be tested for enrichment individually.
#' @return A table reporting enrichment results for each factor given in the "name" column in features.
#' @export
testEnrichment <- function(query, null, features)
{
	# Convert to GRanges
	query.gr <- makeGRanges(query)
	background.gr <- makeGRanges(background)
	features.gr <- makeGRanges(features)

	if(!("name" %in% colnames(values(features.gr))))
	{
		stop("The object given in features must contain a column named \"name\". The set of ranges for each factor type given will be tested separately for enrichment.")
	}

	# Total number of sites in each
	n.query <- length(query.gr)
	n.background <- length(background.gr)

	# Perform master overlap of all queries with all possible features and produce list of counts for each factor type
	message("Computing query to features overlap")
	query.counts <- countFeatures(query.gr, features.gr)

	# Same for the background
	message("Computing query to background overlap")
	background.counts <- countFeatures(background.gr, features.gr)
	names(background.counts) <- c("name","nBackgroundRegions")

	# Combine and fix if there are missings
	#counts <- plyr::join(query.counts, background.counts, by="factor", type="full")
	counts <- merge(query.counts, background.counts, all=T, by="name")
	names(counts) <- c("name", "query", "background")
	counts[is.na(counts[,2]),2] <- 0
	counts[is.na(counts[,3]),3] <- 0

	# Perform enrichment testing
	dotest <- function(count)
	{
		o.query <- count[1]
		o.background <- count[2]
		# Perform binom test
		# Number of sites with feature in query
		x <- o.query
		# Total number of sites in query
		n <- n.query
		# Number of sites with feature in background / Number of sites in background
		p <- o.background / n.background
		if(is.na(p)){P<-0}
		test <- binom.test(x, n, p)
		pv <- test$p.value
		pv
	}
	message("Calculating enrichment test")
	pvs <- apply(counts[,2:3], MARGIN=1, FUN=dotest)

	# Adjust for multiple testing based on the number of features we ran a test for
	pva <- p.adjust(pvs, method="fdr")

	# Output final dataframe
	frac.query <- round(counts$query/n.query,2)
	frac.background <- round(counts$background/n.background,2)
	out <- data.frame(name=counts$name, overlap.query=counts$query, n.query=n.query, overlap.background=counts$background, n.background=n.background, frac.query, frac.background, frac.diff=frac.query-frac.background, p.value=pvs, p.adjusted=pva, sig=pva<0.05)

	out <- out[order(out$frac.diff, decreasing=TRUE),]
	out
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Draw a length-matched pool of sequences from the genome
#'
#' Given a query set of ranges, draw a length-matched pool of sequences. Returned ranges are required to (1) not overlap with each other or the query, (2) not extend off chromosome ends, (3) not extend over assembly gaps as defined in the UCSC "gap" table for the given genome assembly.
#' @param query A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param n Number of times greater than the query set that the size of the returned background pool will be
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @return A GRanges of the background sequences.
#' @export
drawGenomePool <- function(query, n, genome, cachedir, chrs=NULL)
{
	target.gr <- makeGRanges(query)
	if(is.null(chrs))
	{
		chrs <- unique(seqnames(target.gr))
	}

	#target.gr <- dmrs.gr
	lens <- width(target.gr)
	# duplicate the sizes by the number of regions we want to draw per size
	lens <- rep(lens,n)

	# Get chromosome lengths
	chromsizes <- getUCSCTable("chromInfo", genome, cachedir)
	chromsizes <- chromsizes[chromsizes$chrom %in% chrs,]

	# Get gaps GR
	gaps <- getUCSCTable("gap", genome, cachedir)
	gaps.gr <- with(gaps, GRanges(seqnames=chrom, ranges=IRanges(start=chromStart+1, end=chromEnd)))

	# For each size, draw n random non-gap positions

	dodraw <- function(lens)
	{
		chrs <- sample(chrs,length(lens),replace=T)
		dt <- data.table(len=lens,chr=chrs)
		dt <- dt[,list(len=len, start=sample(1:(chromsizes[chromsizes$chrom==chr,]$size),length(len))), by=chr]
		draw.gr <- GRanges(dt$chr,IRanges(start=dt$start,width=dt$len))
		return(draw.gr)
	}
	draw.gr <- dodraw(lens)

	# Exclude the bads and redraw them
	isbad <- TRUE
	iter <- 1
	while(isbad)
	{
		# Check for overlapping regions
		dups <- data.table(as.data.frame(findOverlaps(draw.gr,draw.gr)))
		dups <- dups[queryHits!=subjectHits,]$queryHits
		# Check for those that hit gap regions
		gaps <- data.table(as.data.frame(findOverlaps(draw.gr,gaps.gr)))
		gaps <- gaps[queryHits!=subjectHits,]$queryHits
		# Check for those that extend beyond chr ends
		if(sum(end(draw.gr) > chromsizes[match(as.vector(seqnames(draw.gr)),chromsizes$chrom),]$size)>0){bigs <- seq(1,length(draw.gr))[end(draw.gr) > chromsizes[match(as.vector(seqnames(draw.gr)),chromsizes$chrom),]$size]} else {bigs<-c()}
		#seqlengths(draw.gr) <- chromsizes[match(seqlevels(draw.gr),chromsizes$chrom),]$size
		#1:length(draw.gr)[end(draw.gr) > chromsizes[match(as.vector(seqnames(draw.gr)),chromsizes$chrom),]$size]

		# Need to account for if they are both gaps and dups, otherwise we start building up extra lengths we don't need
		# I fixed this by using unique()
		toget <- unique(c(gaps,dups,bigs))

		# Redraw these widths
		if(length(toget)==0)
		{
			isbad <- FALSE
		} else
		{
			message("Draw iter ",iter," had ", length(dups)," overlapping, ",length(gaps)," in gap regions, and ", length(bigs)," off chromosome ends. Redrawing these...")
			iter <- iter+1
			draw.gr <- c(draw.gr[-c(toget)],dodraw(width(draw.gr[c(toget)])))
		}
	}

	# Make sure the length distribution of the output matches that of the input times the number of draws for each
	if(!all((table(width(target.gr))*n)==table(width(draw.gr)))){stop("Final length tables did not match, the drawing did not work.")}

	draw.gr
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Perform propensity score matching to draw a multi-variate matched set of sequences from a background pool
#'
#' Given a query set and a background pool, draw a set of sequences from the background pool that most closely match the query set with respect to multiple co-variates.
#' @param query A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param pool A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param outdir The function will write a PDF report of matching performance and FASTA files for both the query set and matched null set. Provide a directory in which to save these files.
#' @param formula Formula used for matching. For example: "treat ~ sizeLog + freqCpG". The variable "treat" must always be given as the predicted variable. Combinations of predictors can be selected from the set of: size (length of sequence in bp), sizeLog (log of size - recommended for best matching performance), gc (GC percent), freqCpG (dinucleotide frequency of CpG sites), freqA (frequency of the base A), freqT, freqC, freqG, repeatPer (percent of sequence covered by repeat masked regions), distTSSCenterLogX1 (distance to transcription start sites, log transformed), and distTSECenterLogX1 (distance to transcription end sites, log transformed).
#' @param n Number of times greater than the query the matched null set will be.
#' @param bsg BString genome from which sequence data can be derived. For example, see the "BSgenome.Hsapiens.UCSC.hg19" for the "hg19" genome BSGenome package from Bioconductor. Similar packages exist for other genomes.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @return A list containing the sequences of both the target and pool along with a GRanges of the matched results which can be used as a null set in testEnrichment().
#' @export
doPropMatch <- function(query, pool, outdir=".", formula, n=1, bsg, genome, cachedir)
{
	dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
	target.gr <- makeGRanges(target)
	pool.gr <- makeGRanges(pool)
	formula <- as.formula(formula)

	# Get sequences
	#library(BSgenome.Hsapiens.UCSC.hg19)
	#bsg <- BSgenome.Hsapiens.UCSC.hg19

	# Get sequence/region covariates

	# Return list - seq and meta
	#genome <- "hg19"
	#cachedir <- "~/myucsc"

	message("Meta for seq1")
	seq1 <- goldmine:::getSeqMeta(target.gr, bsg, genome, cachedir)

	# filter out all sequences longer than or shorter than the target from the background pool
	seq2.unfiltered <- pool.gr
	seq2.unfiltered.size <- width(seq2.unfiltered)
	seq2.filt <- seq2.unfiltered[(seq2.unfiltered.size<=max(seq1$meta$size))&(seq2.unfiltered.size>=min(seq1$meta$size))]

	message("Meta for seq2")
	seq2 <- getSeqMeta(seq2.filt, bsg, genome, cachedir)
	bad <- is.na(seq2$meta$gc)
	seq2$meta <- seq2$meta[!bad,]
	seq2$seq <- seq2$seq[!bad]
	
	# set which cols have covars in them
	cols <- 5:ncol(seq1$meta)

	# Do matching
	# Draw Random Starting Order - Use same start order for all runs
	message("Do matching")
	nTotalSeqs <- nrow(seq1$meta) + nrow(seq2$meta)
	index.random <- sample(seq(1:nTotalSeqs),nTotalSeqs, replace=FALSE)
	#formula <- as.formula("treat ~ sizeLog + gc + freqCpG + repeatPer + distTSSCenterLogX1")
	seq.ref <- drawBackgroundSetPropensity(seq1$seq,seq1$meta,seq2$seq,seq2$meta,formula,start.order=index.random,n=n)
	pro.gr <- makeGRanges(seq.ref)

	message("Saving Sequence Sets as FASTA in: ",outdir)
	proseq <- getSeq(bsg, pro.gr)
	names(proseq) <- paste(genome,seqnames(pro.gr),start(pro.gr),end(pro.gr),"+",sep="_")
	names(seq1$seq) <- paste(genome,seqnames(target.gr),start(target.gr),end(target.gr),"+",sep="_")
	#dir.create("output/motif/fasta",showWarnings=FALSE)
	writeXStringSet(proseq,paste0(outdir,"/nullset.fa"))
	writeXStringSet(seq1$seq,paste0(outdir,"/queryset.fa"))

	message("Do plotting")
	# Do plotting
	pdfpath <- paste0(outdir,"/psm.pdf")
	message("Saving Matching Performance PDF: ",pdfpath)
	pdf(file=pdfpath, width=10.5, height=8, paper="USr")
	plotCovarHistogramsOverlap(seq1$meta,seq2$meta,cols,main="Target vs Pool")
	plotCovarHistogramsOverlap(seq1$meta,seq.ref,cols,main="Target vs Matched Reference")
	mymeta <- list(pool=seq2$meta,match=seq.ref)
	plotCovarQQ(seq1$meta, mymeta, cols, plot.ncols=4)
	print(plotCovarDistance(seq1$meta, mymeta, cols))
	dev.off()

	list(query.seq=seq1, null.seq=seq2, ranges.null=pro.gr)
}
# --------------------------------------------------------------------

# ====================================================================

# ====================================================================
# Internal Functions

# --------------------------------------------------------------------
# Given two sets of genomic ranges, count the number of each type of features in the query
countFeatures <- function(query.gr, features.gr)
{
	fo <- findOverlaps(query.gr, features.gr)
	qh <- queryHits(fo)
	sh <- subjectHits(fo)
	query.hits <- data.frame(chr=seqnames(query.gr[qh]), start=start(query.gr[qh]), end=end(query.gr[qh]), name=features.gr[sh]$name)
	query.hits$queryRegion <- with(query.hits,paste(chr, ":", start, "-", end, sep=""))
	#query.counts <- ddply(query.hits, .(factor), summarize, nQueryRegions=length(unique(queryRegion)))
	query.hits <- data.table(query.hits)
	query.counts <- data.frame(query.hits[,list(nQueryRegions=length(unique(queryRegion))),by="name"])
	query.counts
}

# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Make set of distinct color labels for plots
# Input: n=number of colors to output, rand=shuffle color order or not
# Output: Vector of hex color codes
genColors <- function(n, rand=FALSE)
{
	# Brewer's qualitative palette "Set1" only has 9 values
	# Extrapolate from these to create palettes of any size

	pal <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(n)
	if(rand==TRUE){pal <- sample(pal)}

	pal
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Clean-up of ggplot defaults (remove grid)
# Input:
# Output:
ggplot.clean <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = grid::unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Calculate GC content of a DNAString
# Input: DNAStringSet object
# Output: Vector of GC contents for each sequence in the DNAStringSet
getGC <- function(seq)
{
	g <- alphabetFrequency(seq)[,3]
	c <- alphabetFrequency(seq)[,2]
	a <- alphabetFrequency(seq)[,1]
	t <- alphabetFrequency(seq)[,4]

	gc <- (g+c) / (a+t+g+c)
	gc
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Calculate covariates for each sequence in a DNAStringSet
# @param myseq DNAStringSet object of the sequence set
# @return dataframe with standard covariates added
getSeqMeta <- function(ranges,bsgenome,genome,cachedir)
{
	# Input as GRanges
	gr <- makeGRanges(ranges)

	# Get sequence from BSGenome
	seq <- getSeq(bsgenome, gr)
	names(seq) <- 1:length(seq)

	# name, size and gc - easy
	name <- names(seq)
	#name <- 1:length(seq)
	size <- width(seq)
	gc <- goldmine:::getGC(seq)

	# add chr/start/end fields parsed from sequence name
	#name.split <- data.frame(do.call('rbind', strsplit(as.character(name),'-',fixed=TRUE)))
	#names(name.split) <- c("chr","start","end")
	#chr <- name.split$chr
	#start <- as.numeric(as.character(name.split$start))
	#end <- as.numeric(as.character(name.split$end))
	chr <- seqnames(gr)
	start <- start(gr)
	end <- end(gr)

	# Add extra variables from annotation
	seq.ranges <- gr

	repeatPer <- goldmine:::getRepeatPercentFast(seq.ranges,genome,cachedir)
	sizeLog <- log10(size)
	#distTSS <- getDistTSS(seq.ranges,ann)
	#distTSE <- getDistTSE(seq.ranges,ann)
	#distTSSCenter <- getDistTSSCenter(seq.ranges,ann)
	#distTSECenter <- getDistTSECenter(seq.ranges,ann)
	distTSSCenterLogX1 <- log10(goldmine:::getDistTSSCenter(seq.ranges,genome,cachedir)+1)
	distTSECenterLogX1 <- log10(getDistTSECenter(seq.ranges,genome,cachedir)+1)
	freqCpG <- getFreqCpG(seq)

	# single nucleotide frequencies
	mono <- alphabetFrequency(seq)
	freqA <- mono[,1]/size
	freqT <- mono[,4]/size
	freqC <- mono[,2]/size
	freqG <- mono[,3]/size

	# combine into our covariate dataframe
	#seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer, distTSS, distTSSCenter, distTSSCenterLogX1, distTSE, distTSECenter, distTSECenterLogX1)
	seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, freqA, freqT, freqC, freqG, repeatPer, distTSSCenterLogX1, distTSECenterLogX1)
	#seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer)

	# Return list with seq and then df of metadata
	ret <- list(seq=seq,meta=seq.meta)
	ret
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Draw a matched reference set from a reference pool
#
# Use propensity score matching to create a covariate-matched reference set. Note that the matching function is sensitive to the starting order of the input data. This order is required as a variable so it can be fixed between runs.
#
# @param target.seq \code{DNAStringSet} object of the target set
# @param target.meta \code{data.frame} object with the covariates of the target set
# @param pool.seq \code{DNAStringSet} object of the reference pool to draw covariate matched reference set from
# @param pool.meta \code{data.frame} object with the covariates
# @param formula an \code{as.formula} object for the regression used to generate propensity scores
# @param start.order vector of starting order for matching (must be a sequence of integers in any order from 1 to the total number of sequences in both target.seq and pool.seq)
# @return \code{DNAStringSet} object of a covariate-matched reference set
drawBackgroundSetPropensity <- function(target.seq, target.meta, pool.seq, pool.meta, formula, start.order, n=1)
{
	# setting binary value for group assignment
	target.meta$treat <- 1
	pool.meta$treat <- 0
	all.meta <- rbind(target.meta, pool.meta)

	# randomize sort order - order can bias when Match(..., replace=FALSE)
	all.meta.shuffle <- all.meta[start.order,]

	# run logistic model
	lrm.out <- rms::lrm(formula, data=all.meta.shuffle)

	# obtain values
	lrm.out.fitted <- rms:::predict.lrm(lrm.out,type="fitted")

	# if model didn't fit, exclude these from matching
	bad <- is.na(lrm.out.fitted)
	all.meta.shuffle <- all.meta.shuffle[!bad,]
	lrm.out.fitted <- lrm.out.fitted[!bad]

	# match
	rr <- Matching::Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.fitted, M=n, version="standard", replace=FALSE)
	#summary(rr)

	# make new sequence set
	matched.meta <- all.meta.shuffle[rr$index.control,]
	m <- match(as.character(matched.meta$name),names(pool.seq))
	seq.resamp <- pool.seq[m]
	ret <- pool.meta[m,]
	ret
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Perform binomial test of enrichment for motifs using counts of occurrences in two sequence sets
#
# The *.counts matrix objects must first be generated using \code{\link{calcMotifCounts}}. The current implementation only considers if a sequence has at least one occurrence of the motif or not, and does not account for or weight multiple occurrences of a motif in a single sequence. The contingency table is simply based on the number of sequences which contain at least one occurrence of each motif.
#
# @param seq1.counts output object (matrix) from \code{\link{calcMotifCounts}} for first sequence set
# @param seq1.nSeqs number of sequences in first sequence set
# @param seq2.counts output object from \code{\link{calcMotifCounts}} for second sequence set
# @param seq2.nSeqs number of sequences in second sequence set
# @return dataframe of output results including p-values (unadjusted)
calcEnrichmentBinom <- function(seq1.counts,seq1.nSeqs,seq2.counts,seq2.nSeqs)
{
	#input: count matrices (sequences vs motifs) for two runs of FIMO, number of seqs for each set
	#output: enrichment p-values for each motif

	#calculate frequencies for each motif in each set of sequences
	#consider each occurance of motif only once: convert counts >1 to be 1 so we can easily sum each row to get counts
	makeBinary <- function(value)
	{
		if(value>1)
		{
			value <- 1
		}
		value
	}
	counts1.bin <- apply(seq1.counts,MARGIN=c(1,2),FUN=makeBinary)	
	counts2.bin <- apply(seq2.counts,MARGIN=c(1,2),FUN=makeBinary)

	pvalues <- foreach(i=1:nrow(counts1.bin),.combine=rbind) %dopar%
	{
		#Using binom.test:
		#x = num successes = total count of sequences with one or more instance of motif found by FIMO
		myX <- sum(counts1.bin[i,])
		#n = num trials = number of sequences
		myN <- seq1.nSeqs
		#p = prob of success = count of seqs with >=1 instance in background divided by total seqs in background	
		if(rownames(counts1.bin)[1] %in% rownames(counts2.bin))
		{
			#if motif shows up in background set
			#extract row with matching motif name and calculate frequency
			myP <- sum(counts2.bin[rownames(counts2.bin) %in% rownames(counts1.bin)[i],])/seq2.nSeqs
		} else
		{
			#if motif did not show up in background set
			#print("Instance w/o same motif in both")
			myP <- 0
		}

		pValue <- binom.test(x=myX, n=myN, p=myP, alternative="greater")$p.value

		data.frame(motif=rownames(counts1.bin)[i],pvalue=pValue,percent_seqs=((myX/myN)*100))
	}

	pvalues
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Plot histograms in a grid for arbitrary number of variables
#
# Plots non-overlapping single histograms in a grid for all covariate data in a dataframe.
#
# @param seq.meta data.frame of sequence covariates
# @param cols vector of which columns to plot histograms for from seq.meta
# @return plot sent to current graphics device
plotCovarHistograms <- function(seq.meta,cols)
{
	# do we need arguments to take filtering and breaks options?

	# calculate histograms
	hists <- foreach(i=1:length(cols)) %do%
	{
		hist(seq.meta[,cols[i]],plot=FALSE)
	}

	ggplot.hist <- function(h)
	{
		plot.data <- data.frame(bin=h$mids,freq=h$count/sum(h$count))
		ggplot(plot.data, aes(x=bin,y=freq)) + geom_bar(data=plot.data, stat="identity",alpha=0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = grid::unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
	}

	p <- list()
	for(i in 1:length(cols))
	{
  		p[[i]] <- ggplot.hist(hists[[i]]) + labs(title=names(seq.meta)[cols[i]])
	}
	do.call(grid.arrange,p)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Plot overlapping histograms for any number of variables from 2 sets
#
# Plots a grid of overlapping histograms. Data for seq1 will appear in red and seq2 in blue. The region where the distributions will appear in purple.
#
# @param seq1.meta dataframe of covariates from first distribution
# @param seq2.meta dataframe of covariates from second distribution
# @param cols which columns in the *.meta dataframes contain covariate data to plot
# @param plot.ncols number of columns in the plotted grid
# @param main title for the grid of plots (useful if you want to put on the formula you used to generate seq2.meta)
# @return plot to active graphics device
plotCovarHistogramsOverlap <- function(seq1.meta,seq2.meta,cols,plot.ncols=3, main="")
{
	# do we need arguments to take filtering and breaks options?
	name1 <- "Target"
	name2 <- "Reference"

	getbreaks <- function(i)
	{
		scale.max <- max(c(seq1.meta[,cols[i]],seq2.meta[,cols[i]]))
		scale.min <- min(c(seq1.meta[,cols[i]],seq2.meta[,cols[i]]))
		scale.bins <- 20
		seq(from=scale.min,to=scale.max,length.out=scale.bins)
	}
	breaks <- lapply(1:length(cols),getbreaks)

	# calculate histograms
	hist1 <- function(i)
	{
		hist(seq1.meta[,cols[i]],breaks=breaks[[i]],plot=FALSE)
	}
	hists1 <- lapply(1:length(cols),hist1)

	hist2 <- function(i)
	{
		hist(seq2.meta[,cols[i]],breaks=breaks[[i]],plot=FALSE)
	}
	hists2 <- lapply(1:length(cols),hist2)

	ggplot.hist.overlap <- function(h1,h2)
	{
		plot.data <- data.frame(bin=h1$mids,freq=h1$count/sum(h1$count),seq=name1)
		plot.data <- rbind(plot.data,data.frame(bin=h2$mids,freq=h2$count/sum(h2$count),seq=name2))
		ggplot(plot.data, aes(x=bin,y=freq,fill=seq)) + geom_bar(data=subset(plot.data,seq == name1), stat="identity",alpha=0.5) + geom_bar(data=subset(plot.data,seq == name2), stat="identity",alpha=0.5) + scale_fill_manual("", values = c("#007FFF","#FF007F")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = grid::unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
	}

	p1 <- ggplot.hist.overlap(hists1[[1]], hists2[[1]]) + theme(legend.position="bottom")
	tmp <- ggplot_gtable(ggplot_build(p1))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]

	p <- list()
	for(i in 1:length(cols))
	{
  		p[[i]] <- ggplot.hist.overlap(hists1[[i]], hists2[[i]]) + labs(title=names(seq1.meta)[cols[i]]) + theme(legend.position = "none")
	}
	p <- c(p,list(ncol=plot.ncols, main=main))

	g <- do.call(gridExtra::arrangeGrob,p)
	gridExtra::grid.arrange(g, legend, heights=grid::unit(c(7.5,0.5),"in"),nrow=2,ncol=1)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Plot QQ plots in a grid for arbitrary number of variables
#
# Creates QQ plots to compare two distributions for any number of variables
#
# @param orig.meta data frame from the original distribution
# @param list.meta list of data frames for all the distributions to compare to for each variable
# @param cols which columns have covariates to plot from the above dataframes
# @param plot.ncols how many columns the plotted grid should have
# @return plot to active graphics device
plotCovarQQ <- function(orig.meta,list.meta,cols,plot.ncols=3)
{
	ggplot.qq <- function(d,i)
	{
		ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences), size=2, stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background") + geom_abline(slope = 1, intercept=0) + labs(title=names(orig.meta)[cols[i]]) + scale_colour_manual(values = genColors(length(list.meta)))
	}

	getplotdata <- function(i)
	{
		ret <- lapply(1:length(list.meta),function(u) data.frame(as.data.frame(qqplot(orig.meta[,cols[i]], list.meta[[u]][,cols[i]], plot.it=FALSE)),Sequences=names(list.meta)[u]))
		ret <- do.call("rbind",ret)
		ret
	}
	plot.data <- lapply(1:length(cols),getplotdata)

	p1 <- ggplot.qq(plot.data[[1]]) + theme(legend.position="right")
	tmp <- ggplot_gtable(ggplot_build(p1))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]

	p <- list()
	for(i in 1:length(cols))
	{
		d <- plot.data[[i]]
  		p[[i]] <- ggplot.qq(d,i) + theme(legend.position="none")
	}
	p <- c(p,list(ncol=plot.ncols))

	g <- do.call(gridExtra::arrangeGrob,p)
	#grid.arrange(g, legend, widths=unit(c(7.5,0.5),"in"), main="QQ Plots",nrow=1,ncol=2)
	gridExtra::grid.arrange(g, main="QQ Plots")
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Plot horizontal graph of covariate distance from different propensity models
#
# Horizontal graph plots a point for each variable that represents the distance between that variable's value in orig.meta and each of the dataframes in list.meta
#
# @param orig.meta dataframe of covariates from the target set you are trying to match
# @param list.meta a list of dataframes of covariates from other sets you want to compare to the target set
# @param cols vector of which columns to use from the dataframes above
# @return plots to active graphics device
plotCovarDistance <- function(orig.meta,list.meta,cols)
{
	stddist <- function(d1, d2)
	{
		#(mean(d1)-mean(d2))/(sd(c(d1,d2))/2)
		(100*abs(mean(d1)-mean(d2)))/sqrt(((sd(d1)^2)+(sd(d2)^2))/2)
	}

	# calculate distances for each variable
	getdists <- function(i)
	{
		ret <- sapply(1:length(cols),function(u) stddist(orig.meta[,cols[u]],list.meta[[i]][,cols[u]]))
		ret
	}
	dists <- lapply(1:length(list.meta),getdists)
	dists <- do.call("rbind", dists)

	colnames(dists) <- names(orig.meta)[cols]
	rownames(dists) <- names(list.meta)

	plot.data <- reshape::melt.matrix(dists)
	names(plot.data) <- c("matching","variable","stddist")

	plot.data.sub <- plot.data[plot.data$matching=="pool",]
	plot.data.sub$variable <- as.character(plot.data.sub$variable)
	mylevs <- plot.data.sub[sort(plot.data.sub$stddist, index.return=TRUE)$ix,]$variable

	#mylevs <- levels(reorder(x=plot.data[plot.data$matching=="pool",]$variable,X=plot.data[plot.data$matching=="pool",]$stddist, order=FALSE))

	plot.data$matching <- factor(plot.data$matching,levels=names(list.meta))

	plot.data$variable <- factor(plot.data$variable,levels=mylevs)

	ggplot(plot.data, aes(x=variable,y=stddist,col=matching)) + geom_point(data=plot.data,size=3) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linetype=3, colour="grey50"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = grid::unit(0.8, "lines"), axis.line = element_line(colour = "grey50"), axis.text=element_text(colour="black")) + geom_abline(intercept=0,slope=0,col="grey50") + coord_flip() + labs(main="Covariate Balance",y="Standardized Distance") + scale_colour_manual(values = genColors(length(list.meta)))
}
# --------------------------------------------------------------------

# ====================================================================

