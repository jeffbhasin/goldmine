# Functions for testing enrichment between sets of genomic ranges and generating background sets

# -----------------------------------------------------------------------------
#' Make a GRanges from a data.frame or data.table with the fields "chr", "start", and "end"
#'
#' Given a data.frame or data.table with the columns "chr", "start", and "end", a GenomicRanges (GRanges) object will be created. All other columns will be passed on as metadata. If the input is already a GRanges, it is simply returned. If the column "strand" exists, it will be set as the strand.
#' @param obj A data.frame or data.table with columns "chr", "start", and "end" and any other columns
#' @param strand Use the information in the "strand" column to set strand in the GRanges, if it is present.
#' @return A GRanges made from the data in obj.
#' @export
makeGRanges <- function(obj, strand=F)
{
	if(class(obj)[1]=="GRanges")
	{
		# Return if it is already a GRanges
		obj
	} else if(is.data.frame(obj)==TRUE)
	{
		# Make GRanges if it is a data.frame
		if(sum(c("chr", "start", "end") %in% colnames(obj))!=3)
		{
			stop("Could not find columns named \"chr\", \"start\", and \"end\" in input data.frame")
		}
		ret <- with(obj,GRanges(seqnames=chr,IRanges(start,end)))
		if(("strand" %in% colnames(obj))&(strand==T))
		{
			strand(ret) <- obj$strand
		}
		skipcols <- c("chr","start","end","strand","width","element")
		if(class(obj)[1]=="data.table")
		{
			values(ret) <- obj[,!(colnames(obj) %in% skipcols),with=F]
		} else
		{
			values(ret) <- obj[,!(colnames(obj) %in% skipcols)]
		}
		ret
	} else
	{
		# Some bad input, throw error
		stop("GRanges, data.frame, or data.table object required as input")
	}
}
# -----------------------------------------------------------------------------

# Given two sets of genomic ranges, count the number of each type of features in the query
countFeatures <- function(query.gr, features.gr)
{
	fo <- findOverlaps(query.gr, features.gr)
	qh <- queryHits(fo)
	sh <- subjectHits(fo)
	query.hits <- data.frame(chr=seqnames(query.gr[qh]), start=start(query.gr[qh]), end=end(query.gr[qh]), factor=features.gr[sh]$factor)
	query.hits$queryRegion <- with(query.hits,paste(chr, ":", start, "-", end, sep=""))
	#query.counts <- ddply(query.hits, .(factor), summarize, nQueryRegions=length(unique(queryRegion)))
	query.hits <- data.table(query.hits)
	query.counts <- data.frame(query.hits[,list(nQueryRegions=length(unique(queryRegion))),by="factor"])
	query.counts
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Make a GRanges from a data.frame or data.table with the fields "chr", "start", and "end"
#'
#' Given a data.frame or data.table with the columns "chr", "start", and "end", a GenomicRanges (GRanges) object will be created. All other columns will be passed on as metadata. If the input is already a GRanges, it is simply returned.
#' @param obj A data.frame or data.table with columns "chr", "start", and "end" and any other columns
#' @return A GRanges made from the data in obj.
#' @export
testEnrichment <- function(query, background, features)
{
	# Convert to GRanges
	query.gr <- makeGRanges(query)
	background.gr <- makeGRanges(background)
	features.gr <- makeGRanges(features)

	if(!("factor" %in% colnames(values(features.gr))))
	{
		stop("The object given in features must contain a column named \"factor\". The set of ranges for each factor type given will be tested separately for enrichment.")
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
	names(background.counts) <- c("factor","nBackgroundRegions")

	# Combine and fix if there are missings
	#counts <- plyr::join(query.counts, background.counts, by="factor", type="full")
	counts <- merge(query.counts, background.counts, all=T, by="factor")
	names(counts) <- c("factor", "query", "background")
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
	out <- data.frame(factor=counts$factor, overlap.query=counts$query, n.query=n.query, overlap.background=counts$background, n.background=n.background, frac.query, frac.background, frac.diff=frac.query-frac.background, p.value=pvs, p.adjusted=pva, sig=pva<0.05)

	out <- out[order(out$frac.diff, decreasing=TRUE),]
	out
}
# -----------------------------------------------------------------------------

# Draw random regions from anywhere in the genome with the same length distribution as the target set
# Drawn output will not be in assembly gaps, will not run off chromosome ends, and will not overlap with each other
drawBackgroundSetGenomic <- function(n, target.gr, genome, cachedir, chrs)
{
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

	# Make empty draws.gr
	draws.gr <- GRanges()

	# For each size, draw n random non-gap positions
	for(len in lens)
	{
		message("Finding length ", len)

		isbad <- TRUE
		while(isbad)
		{
			# Draw a random chr
			chr <- sample(chrs,1)
	
			# Make bads GR of all regions within the chromosome size we don't want to draw from (gaps and previous drawn regions)
			bad.gr <- reduce(c(gaps.gr, draws.gr))

			# Draw random start position, up to the closest we can get to the chr end without running off based on our size
			rand.start <- sample(1:(chromsizes[chromsizes$chrom==chr,]$size-len),1)
			rand.gr <- GRanges(seqnames=chr, ranges=IRanges(start=rand.start, width=len))
	
			# Test if we drew into a bad region - need to redraw if we did
			isbad <- rand.gr %over% bad.gr
			message("Draw bad? ",isbad)
		}
		draws.gr <- suppressWarnings(c(draws.gr, rand.gr))

	}
	
	draws.gr
}

# Do checks at the end so we can do draws in parallel
drawBackgroundSetGenomicFast <- function(n, target.gr, genome, cachedir, chrs)
{
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