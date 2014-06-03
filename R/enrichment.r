# Functions for testing enrichment between sets of genomic ranges and generating background sets

# -----------------------------------------------------------------------------
#' Make a GRanges from a data.frame or data.table with the fields "chr", "start", and "end"
#'
#' Given a data.frame or data.table with the columns "chr", "start", and "end", a GenomicRanges (GRanges) object will be created. All other columns will be passed on as metadata. If the input is already a GRanges, it is simply returned.
#' @param obj A data.frame or data.table with columns "chr", "start", and "end" and any other columns
#' @return A GRanges made from the data in obj.
#' @export
makeGRanges <- function(obj)
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
	counts <- join(query.counts, background.counts, by="factor", type="full")
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
