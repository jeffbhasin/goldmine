# #############################################################################
# Finds differentially enriched motifs
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013
# #############################################################################

# =============================================================================
# Packages and Globals
library(ShortRead)
library(plyr)
library(reshape)
library(foreach)
library(stringr)
library(ggplot2)
library(gridExtra)
library(Matching) # for Match()
library(rms)
library(doMC)
library(gplots) # for heatmap.2()
library(VennDiagram)
# =============================================================================

# =============================================================================
# Utility

# -----------------------------------------------------------------------------
# Make set of distinct color labels for plots
# Input: n=number of colors to output, rand=shuffle color order or not
# Output: Vector of hex color codes
genColors <- function(n, rand=FALSE)
{
	# Brewer's qualitative palette "Set1" only has 9 values
	# Extrapolate from these to create palettes of any size

	pal <- colorRampPalette(brewer.pal(9,"Set1"))(n)
	if(rand==TRUE){pal <- sample(pal)}

	pal
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Clean-up of ggplot defaults (remove grid)
# Input:
# Output:
ggplot.clean <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# =============================================================================


# =============================================================================
# MEME Suite Interaction

# -----------------------------------------------------------------------------
#' Wrapper to call FIMO using system()
#'
#' If FIMO from the MEME Suite is installed and in the current PATH, this provides an easy interface
#' to run it from R.
#'
#' @param out.path path where output will be written
#' @param fasta.path path to FASTA file of input sequences
#' @param motifs.path path to MEME format file containing motif database to use
#' @export
runFIMO <- function(out.path,fasta.path,motifs.path)
{
	system(paste("fimo -oc ",out.path," ",motifs.path," ",fasta.path,sep=""))
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Read in a FIMO output file
#'
#' Wrapper of read.table() with correct options for reading in a FIMO output text file.
#'
#' @param fimo.out.path path to FIMO output text file
#' @export
readFIMO <- function(fimo.out.path)
{
	read.table(file=fimo.out.path,header=TRUE,comment.char="",sep="\t")
}
# -----------------------------------------------------------------------------
# =============================================================================

# =============================================================================
# Background Sequence Generation

#' Calculate covariates for each sequence in a DNAStringSet
#' @param myseq DNAStringSet object of the sequence set
#' @return dataframe with standard covariates added
#' @export
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
	gc <- getGC(seq)

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

	repeatPer <- getRepeatPercentFast(seq.ranges,genome,cachedir)
	sizeLog <- log10(size)
	#distTSS <- getDistTSS(seq.ranges,ann)
	#distTSE <- getDistTSE(seq.ranges,ann)
	#distTSSCenter <- getDistTSSCenter(seq.ranges,ann)
	#distTSECenter <- getDistTSECenter(seq.ranges,ann)
	distTSSCenterLogX1 <- log10(getDistTSSCenter(seq.ranges,ann)+1)
	distTSECenterLogX1 <- log10(getDistTSECenter(seq.ranges,ann)+1)
	freqCpG <- getFreqCpG(seq)

	# combine into our covariate dataframe
	#seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer, distTSS, distTSSCenter, distTSSCenterLogX1, distTSE, distTSECenter, distTSECenterLogX1)
	seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer, distTSSCenterLogX1, distTSECenterLogX1)
	#seq.meta <- data.frame(name, chr, start, end, size, sizeLog, gc, freqCpG, repeatPer)

	# Return list with seq and then df of metadata
	ret <- list(seq=seq,meta=seq.meta)
	ret
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Draw a matched reference set from a reference pool
#'
#' Use propensity score matching to create a covariate-matched reference set. Note that the matching function is sensitive to the starting order of the input data. This order is required as a variable so it can be fixed between runs.
#'
#' @param target.seq \code{DNAStringSet} object of the target set
#' @param target.meta \code{data.frame} object with the covariates of the target set
#' @param pool.seq \code{DNAStringSet} object of the reference pool to draw covariate matched reference set from
#' @param pool.meta \code{data.frame} object with the covariates
#' @param formula an \code{as.formula} object for the regression used to generate propensity scores
#' @param start.order vector of starting order for matching (must be a sequence of integers in any order from 1 to the total number of sequences in both target.seq and pool.seq)
#' @return \code{DNAStringSet} object of a covariate-matched reference set
#' @export
drawBackgroundSetPropensity <- function(target.seq, target.meta, pool.seq, pool.meta, formula, start.order)
{
	# setting binary value for group assignment
	target.meta$treat <- 1
	pool.meta$treat <- 0
	all.meta <- rbind(target.meta, pool.meta)

	# randomize sort order - order can bias when Match(..., replace=FALSE)
	all.meta.shuffle <- all.meta[start.order,]

	# run logistic model
	lrm.out <- lrm(formula, data=all.meta.shuffle)

	# obtain values
	lrm.out.fitted <- rms:::predict.lrm(lrm.out,type="fitted")

	# match
	rr <- Match(Y=NULL, Tr=all.meta.shuffle$treat, X=lrm.out.fitted, M=1, version="standard", replace=FALSE)
	#summary(rr)

	# make new sequence set
	matched.meta <- all.meta.shuffle[rr$index.control,]
	m <- match(as.character(matched.meta$name),names(pool.seq))
	seq.resamp <- pool.seq[m]
	ret <- pool.meta[m,]
	ret
}
# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Plots

# -----------------------------------------------------------------------------
#' Plot histograms in a grid for arbitrary number of variables
#'
#' Plots non-overlapping single histograms in a grid for all covariate data in a dataframe.
#'
#' @param seq.meta data.frame of sequence covariates
#' @param cols vector of which columns to plot histograms for from seq.meta
#' @return plot sent to current graphics device
#' @export
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
		ggplot(plot.data, aes(x=bin,y=freq)) + geom_bar(data=plot.data, stat="identity",alpha=0.8) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
	}

	p <- list()
	for(i in 1:length(cols))
	{
  		p[[i]] <- ggplot.hist(hists[[i]]) + labs(title=names(seq.meta)[cols[i]])
	}
	do.call(grid.arrange,p)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Plot overlapping histograms for any number of variables from 2 sets
#'
#' Plots a grid of overlapping histograms. Data for seq1 will appear in red and seq2 in blue. The region where the distributions will appear in purple.
#'
#' @param seq1.meta dataframe of covariates from first distribution
#' @param seq2.meta dataframe of covariates from second distribution
#' @param cols which columns in the *.meta dataframes contain covariate data to plot
#' @param plot.ncols number of columns in the plotted grid
#' @param main title for the grid of plots (useful if you want to put on the formula you used to generate seq2.meta)
#' @return plot to active graphics device
#' @export
plotCovarHistogramsOverlap <- function(seq1.meta,seq2.meta,cols,plot.ncols=3, main="")
{
	# do we need arguments to take filtering and breaks options?
	name1 <- "Target"
	name2 <- "Reference"

	breaks <- foreach(i=1:length(cols)) %do%
	{
		scale.max <- max(c(seq1.meta[,cols[i]],seq2.meta[,cols[i]]))
		scale.min <- min(c(seq1.meta[,cols[i]],seq2.meta[,cols[i]]))
		scale.bins <- 20
		seq(from=scale.min,to=scale.max,length.out=scale.bins)
	}

	# calculate histograms
	hists1 <- foreach(i=1:length(cols)) %do%
	{
		hist(seq1.meta[,cols[i]],breaks=breaks[[i]],plot=FALSE)
	}
	hists2 <- foreach(i=1:length(cols)) %do%
	{
		hist(seq2.meta[,cols[i]],breaks=breaks[[i]],plot=FALSE)
	}

	ggplot.hist.overlap <- function(h1,h2)
	{
		plot.data <- data.frame(bin=h1$mids,freq=h1$count/sum(h1$count),seq=name1)
		plot.data <- rbind(plot.data,data.frame(bin=h2$mids,freq=h2$count/sum(h2$count),seq=name2))
		ggplot(plot.data, aes(x=bin,y=freq,fill=seq)) + geom_bar(data=subset(plot.data,seq == name1), stat="identity",alpha=0.5) + geom_bar(data=subset(plot.data,seq == name2), stat="identity",alpha=0.5) + scale_fill_manual("", values = c("#007FFF","#FF007F")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"))
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

	g <- do.call(arrangeGrob,p)
	grid.arrange(g, legend, heights=unit(c(7.5,0.5),"in"),nrow=2,ncol=1)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Plot QQ plots in a grid for arbitrary number of variables
#'
#' Creates QQ plots to compare two distributions for any number of variables
#'
#' @param orig.meta data frame from the original distribution
#' @param list.meta list of data frames for all the distributions to compare to for each variable
#' @param cols which columns have covariates to plot from the above dataframes
#' @param plot.ncols how many columns the plotted grid should have
#' @return plot to active graphics device
#' @export
plotCovarQQ <- function(orig.meta,list.meta,cols,plot.ncols=3)
{
	ggplot.qq <- function(d)
	{
		ggplot(d) + geom_point(aes(x=x, y=y,color=Sequences), size=2, stat = "identity", position = "identity", ) + ggplot.clean() + labs(x="Target", y="Background") + geom_abline(slope = 1, intercept=0) + labs(title=names(orig.meta)[cols[i]]) + scale_colour_manual(values = genColors(length(list.meta)))
	}

	plot.data <- foreach(i=1:length(cols)) %do%
	{
		foreach(u=1:length(list.meta),.combine="rbind") %do%
		{
			data.frame(as.data.frame(qqplot(orig.meta[,cols[i]], list.meta[[u]][,cols[i]], plot.it=FALSE)),Sequences=names(list.meta)[u])
		}
	}

	p1 <- ggplot.qq(plot.data[[1]]) + theme(legend.position="right")
	tmp <- ggplot_gtable(ggplot_build(p1))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]

	p <- list()
	for(i in 1:length(cols))
	{
		d <- plot.data[[i]]
  		p[[i]] <- ggplot.qq(d) + theme(legend.position="none")
	}
	p <- c(p,list(ncol=plot.ncols))

	g <- do.call(arrangeGrob,p)
	#grid.arrange(g, legend, widths=unit(c(7.5,0.5),"in"), main="QQ Plots",nrow=1,ncol=2)
	grid.arrange(g, main="QQ Plots")
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Plot horizontal graph of covariate distance from different propensity models
#'
#' Horizontal graph plots a point for each variable that represents the distance between that variable's value in orig.meta and each of the dataframes in list.meta
#'
#' @param orig.meta dataframe of covariates from the target set you are trying to match
#' @param list.meta a list of dataframes of covariates from other sets you want to compare to the target set
#' @param cols vector of which columns to use from the dataframes above
#' @return plots to active graphics device
#' @export
plotCovarDistance <- function(orig.meta,list.meta,cols)
{
	stddist <- function(d1, d2)
	{
		#(mean(d1)-mean(d2))/(sd(c(d1,d2))/2)
		(100*abs(mean(d1)-mean(d2)))/sqrt(((sd(d1)^2)+(sd(d2)^2))/2)
	}

	# calculate distances for each variable
	dists <- foreach(i=1:length(list.meta),.combine="rbind") %do%
	{
		vec <- foreach(u=1:length(cols),.combine="c") %do%
		{		
			stddist(orig.meta[,cols[u]],list.meta[[i]][,cols[u]])
		}
	}
	colnames(dists) <- names(orig.meta)[cols]
	rownames(dists) <- names(list.meta)

	plot.data <- melt(dists)
	names(plot.data) <- c("matching","variable","stddist")

	plot.data.sub <- plot.data[plot.data$matching=="pool",]
	plot.data.sub$variable <- as.character(plot.data.sub$variable)
	mylevs <- plot.data.sub[sort(plot.data.sub$stddist, index.return=TRUE)$ix,]$variable

	#mylevs <- levels(reorder(x=plot.data[plot.data$matching=="pool",]$variable,X=plot.data[plot.data$matching=="pool",]$stddist, order=FALSE))

	plot.data$matching <- factor(plot.data$matching,levels=names(list.meta))

	plot.data$variable <- factor(plot.data$variable,levels=mylevs)

	ggplot(plot.data, aes(x=variable,y=stddist,col=matching)) + geom_point(data=plot.data,size=3) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linetype=3, colour="grey50"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.key.size = unit(0.8, "lines"), axis.line = element_line(colour = "grey50"), axis.text=element_text(colour="black")) + geom_abline(intercept=0,slope=0,col="grey50") + coord_flip() + labs(main="Covariate Balance",y="Standardized Distance") + scale_colour_manual(values = genColors(length(list.meta)))
}
# -----------------------------------------------------------------------------

# =============================================================================

# =============================================================================
# Enrichment Testing

# -----------------------------------------------------------------------------
#' Perform binomial test of enrichment for motifs using counts of occurrences in two sequence sets
#'
#' The *.counts matrix objects must first be generated using \code{\link{calcMotifCounts}}. The current implementation only considers if a sequence has at least one occurrence of the motif or not, and does not account for or weight multiple occurrences of a motif in a single sequence. The contingency table is simply based on the number of sequences which contain at least one occurrence of each motif.
#'
#' @param seq1.counts output object (matrix) from \code{\link{calcMotifCounts}} for first sequence set
#' @param seq1.nSeqs number of sequences in first sequence set
#' @param seq2.counts output object from \code{\link{calcMotifCounts}} for second sequence set
#' @param seq2.nSeqs number of sequences in second sequence set
#' @return dataframe of output results including p-values (unadjusted)
#' @export
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
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Generate counts of motif occurrences in each sequence
#'
#' @param fimo.out dataframe from \code\link{readFIMO}}} dataframe object
#' @param q.cutoff only count a motif if the q-value is less than this cutoff
#'
#' @return matrix with pairwise counts of each motif and each sequence
#' @export
calcMotifCounts <- function(fimo.out, q.cutoff)
{
	#input: dataframe of fimo's text output, q value cutoff
	#output: matrix of frequencies with a row for every motif and column for every sequence

	fimo.out.filtered <- fimo.out[fimo.out$q.value<q.cutoff,]

	fimo.out.filtered$X.pattern.name <- as.character(fimo.out.filtered$X.pattern.name)
	fimo.out.filtered$sequence.name <- as.character(fimo.out.filtered$sequence.name)

	tm <- table(fimo.out.filtered$X.pattern.name,fimo.out.filtered$sequence.name)
}
# -----------------------------------------------------------------------------

# =============================================================================
