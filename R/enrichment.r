# Functions for testing enrichment between sets of genomic ranges and generating background sets

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Draw a length-matched pool of sequences from the genome
#'
#' Given a query set of ranges, draw a length-matched pool of sequences. Returned ranges are required to (1) not overlap with each other or the query, (2) not extend off chromosome ends, (3) not extend over assembly gaps as defined in the UCSC "gap" table for the given genome assembly.
#' @param query A data.frame or data.table with columns "chr", "start", and "end" and any other columns. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based.
#' @param n Number of times greater than the query set that the size of the returned background pool will be
#' @param chrs Vector of chromosome names to draw from. Useful for restricting to canonical chromosomes only, i.e. chr1 to chr22, chrX, and chrY for hg19. If not given, will restrict to the chromosome names present in query.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @return A GRanges of the background sequences.
#' @export
drawGenomePool <- function(query, n, chrs=NULL, genome, cachedir, sync=TRUE)
{
	target.gr <- makeGRanges(query)
	if(is.null(chrs))
	{
		message("Drawing from chrs: ",toString(seqlevels(target.gr))," because these were present in query. Please provide the \"chrs\" option if you wish to draw from additional chromosomes.")
		chrs <- unique(seqnames(target.gr))
	}

	#target.gr <- dmrs.gr
	lens <- width(target.gr)
	# duplicate the sizes by the number of regions we want to draw per size
	lens <- rep(lens,n)

	# Get chromosome lengths
	chromsizes <- getUCSCTable("chromInfo", genome, cachedir, sync=sync)
	chromsizes <- chromsizes[chromsizes$chrom %in% chrs,]

	# Get gaps GR
	gaps <- getUCSCTable("gap", genome, cachedir, sync=sync)
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
