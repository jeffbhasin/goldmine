# #############################################################################
# Annotation, integrated with UCSC genome browser table fetch functions
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013-08-21
# #############################################################################

# =============================================================================
# User-Facing Functions

# -----------------------------------------------------------------------------
#' Annotate regions in mutually exclusive genic, flanking, and intergenic categories
#'
#' Given a set of query regions in a data.frame containing the columns "chr", "start", and "end", returns a new data.frame with the same number of rows as the query and adds annotation columns for the percent of each query region that is genic, flanking, or intergenic, in addition to listing the gene symbols involved for each category. A "call" column is generated which allows plotting a pie chart or stacked bar chart.
#' @param query A data.frame of regions to annotate. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. All additional columns will be retained in the output object.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly.
#' @param flank.bp The size in base pairs for the flanking category. Note that portions of flanks that overlap with gene bodies are no longer considered flanks in this analysis.
#' @param gene.names Include gene names (symbols) in output. Disabling improves speed for large sets of regions (100,000+).
#' @return A data.frame with all rows and columns of the query data.frame with the annotation columns added.
#' @export
annotateSimple <- function(query, genome, cachedir=NULL, flank.bp=1000, gene.names=TRUE)
{
	# Check that query has the needed columns
	if(sum(c("chr", "start", "end") %in% colnames(query))!=3)
	{
		stop("Could not find columns named \"chr\", \"start\", and \"end\" in query data.frame")
	}

	print("Generating annotation regions from UCSC knownGene table")

	# Load the UCSC tables we will need
	knownGene <- getUCSCTable("knownGene", genome, cachedir)
	kgXref <- getUCSCTable("kgXref", genome, cachedir)
	knownCanonical <- getUCSCTable("knownCanonical", genome, cachedir)
	chromInfo <- getUCSCTable("chromInfo", genome, cachedir)

	# join in gene symbols
	names(kgXref)[1] <- "name"
	kg <- join(knownGene, kgXref, by="name", type="left", match="first")

	# join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	kg <- join(kg,kgCanon.sub, by="name", type="left", match="first")
	kg[is.na(kg$canonical),]$canonical <- 0

	# add row number field for future rejoins
	kg$srow <- 1:nrow(kg)

	# add row number field to input for future rejoins
	query$qrow <- 1:nrow(query)

	# Turn input data into GRanges
	input.gr <- GRanges(seqnames=query$chr, ranges=IRanges(start=query$start, end=query$end), qrow=query$qrow)
	seqlengths(input.gr) <- chromInfo[match(names(seqlengths(input.gr)), chromInfo$chr),]$size

	# Make genic GRanges
	genic.gr <- with(kg, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1,end=txEnd), srow=srow))
	genic.full.gr <- genic.gr
	genic.gr <- reduce(genic.gr)
	seqlengths(genic.gr) <- chromInfo[match(names(seqlengths(genic.gr)), chromInfo$chr),]$size

	# Make flank GRanges
	flank.gr <- suppressWarnings(c(flank(genic.gr, flank.bp, start=T), flank(genic.gr, flank.bp, start=F)))
	flank.full.gr <- suppressWarnings(c(flank(genic.full.gr, flank.bp, start=T), flank(genic.full.gr, flank.bp, start=F)))
	flank.gr <- reduce(flank.gr)
	flank.gr <- setdiff(flank.gr, genic.gr)

	# Make intergenic GRanges
	intergenic.gr <- GRanges(seqnames=chromInfo$chrom, ranges=IRanges(start=1, end=chromInfo$size))
	seqlengths(intergenic.gr) <- chromInfo[match(names(seqlengths(intergenic.gr)), chromInfo$chr),]$size
	intergenic.gr <- setdiff(intergenic.gr, genic.gr)
	intergenic.gr <- setdiff(intergenic.gr, flank.gr)

	# Debug: write BED files to make sure we have the sets correct
	#writeBEDFromGRanges(flank.gr, file="flank.bed")
	#writeBEDFromGRanges(genic.gr, file="genic.bed")
	#writeBEDFromGRanges(intergenic.gr, file="intergenic.bed")

	# Debug: Check our set math, should all return empty sets (mutually exclusive ranges)
	#intersect(genic.gr, intergenic.gr)
	#intersect(genic.gr, flank.gr)
	#intersect(flank.gr, intergenic.gr)

	# Annotate with % overlaps and gene names
	print("Computing % Overlaps")
	ann <- data.frame(genic.per=calcPercentOverlap(input.gr, genic.gr), flank.per=calcPercentOverlap(input.gr, flank.gr), intergenic.per=calcPercentOverlap(input.gr, intergenic.gr))

	# Internal function to return list of gene names (no duplicates that are in the reported genic overlap or belong to all or part of the reported flanking overlap)
	# Will return name number of rows as subject.gr, with blanks reported if no genes overlapped
	# Needs a kg knownGene table so we can pull the names
	getGeneList <- function(query.gr, subject.gr, kg)
	{
		# Perform overlaps
		overlaps.fo <- findOverlaps(query.gr, subject.gr)

		# Make DF with gene symbol for each overlap
		genes.df <- data.frame(qrow=queryHits(overlaps.fo), srow=subjectHits(overlaps.fo))
		genes.df$gene <- kg[genes.df$srow,]$geneSymbol

		# Condense into comma-separated lists, removing duplicate symbols
		genes.dt <- data.table(genes.df)
		setkey(genes.dt, qrow, srow)

		i<-0
		commaGenes <- function(x)
		{
			i<<-i+1
			if((i %% 10000)==0)
			{
				print(paste("Done with genic gene names for ",i, " regions",sep=""))
			}
			#genes <- unique(kg[x,]$geneSymbol)
			genes <- unique(x)
			genes <- paste(genes, collapse=", ")
			genes
		}
		g <- genes.dt[,commaGenes(gene), by=qrow]
		setkey(g, qrow)

		# Add blanks for intergenics
		out <- data.table(qrow=query.gr$qrow)
		setkey(out, qrow)
		out <- merge(out, g, by="qrow", all.x=TRUE, all.y=FALSE)
		if(sum(is.na(out$V1)) > 0)
		{
			out[is.na(V1),]$V1 <- ""
		}
		# Return vector of strings, rows matching input.gr
		out$V1
	}
	# maps IDs to use the flank list
	# to get these, we need both the set reduced/diffed version and the non. We'll intersect with the reduced version then intersect back to the non-reduced version to pull out what genes were involved in that flank.
	getGeneListForFlanks <- function(query.gr, flank.gr, flank.full.gr, kg)
	{
		# Perform overlaps against reduced list
		overlaps.fo <- findOverlaps(query.gr, flank.gr)
		o.f.gr <- flank.gr[subjectHits(overlaps.fo)]
		o.f.gr$qrow <- queryHits(overlaps.fo)

		# Now we have pairs of reduced regions and query regions
		# Want to know which unreduced srow contributed to the reduced version
		# Intersect these regions back to the original to find out!
		o2.fo <- findOverlaps(o.f.gr, flank.full.gr)

		# Make DF with gene symbol for each overlap
		genes.df <- data.frame(qrow=o.f.gr[queryHits(o2.fo)]$qrow, srow=subjectHits(o2.fo))
		genes.df$gene <- kg[flank.full.gr[genes.df$srow]$srow,]$geneSymbol

		# Condense into comma-separated lists, removing duplicate symbols
		genes.dt <- data.table(genes.df)
		setkey(genes.dt, qrow, srow)

		i<-0
		commaGenes <- function(x)
		{
			i<<-i+1
			if((i %% 10000)==0)
			{
				print(paste("Done with flanking gene names for ",i, " regions",sep=""))
			}
			#genes <- unique(kg[flank.full.gr[x,]$srow,]$geneSymbol)
			genes <- unique(x)
			genes <- paste(genes, collapse=", ")
			#print(genes)
			genes
		}
		g <- genes.dt[,commaGenes(gene), by=qrow]
		setkey(g, qrow)

		# Add blanks for intergenics
		out <- data.table(qrow=query.gr$qrow)
		setkey(out, qrow)
		out <- merge(out, g, by="qrow", all.x=TRUE, all.y=FALSE)
		if(sum(is.na(out$V1)) > 0)
		{
			out[is.na(V1),]$V1 <- ""
		}
		# Return vector of strings, rows matching input.gr
		out$V1
	}
	# find nearest gene, return DF of nearest and distance to
	getNearestGene <- function(query.gr, subject.gr, kg)
	{
		near.fo <- as.data.frame(suppressWarnings(nearest(query.gr, subject.gr, select="all")))
		near.fo$dist <- suppressWarnings(distance(query.gr[near.fo$queryHits], subject.gr[near.fo$subjectHits]))

		# Make DF with gene symbol for each overlap
		genes.df <- data.frame(qrow=near.fo$queryHits, srow=near.fo$subjectHits, dist=near.fo$dist)
		genes.df$gene <- kg[subject.gr[genes.df$srow]$srow,]$geneSymbol

		# Condense into comma-separated lists, removing duplicate symbols
		genes.dt <- data.table(genes.df)
		setkey(genes.dt, qrow, srow)

		i<-0
		commaGenes <- function(x)
		{
			i<<-i+1
			if((i %% 10000)==0)
			{
				print(paste("Done with nearest gene names for ",i, " regions",sep=""))
			}
			#genes <- unique(kg[flank.full.gr[x,]$srow,]$geneSymbol)
			genes <- unique(x)
			genes <- paste(genes, collapse=", ")
			#print(genes)
			genes
		}
		g <- genes.dt[,commaGenes(gene), by=qrow]
		setkey(g, qrow)

		# Add blanks for intergenics
		out <- data.table(qrow=query.gr$qrow)
		out <- merge(out, g, by="qrow", all.x=FALSE, all.y=FALSE)
		if(sum(is.na(out$V1)) > 0)
		{
			out[is.na(V1),]$V1 <- ""
		}
		out <- data.frame(out)
		# Add distance in bp
		ddf <- data.frame(qrow=near.fo$queryHits, nearest.bp=near.fo$dist)
		out <- join(out, ddf, by="qrow", type="left", match="first")
		colnames(out) <- c("qrow", "nearest.gene", "nearest.bp")
		out$qrow <- NULL
		# Return df
		out

	}

	# call categories
	print("Calling categories")
	ann$gen <- NA
	ann[ann$genic.per>0,]$gen <- "genic"
	ann$cis <- NA
	ann[ann$flank.per>0,]$cis <- "flank"
	ann$int <- NA
	ann[ann$intergenic.per>0,]$int <- "intergenic"

	ann$call <- paste(ann$gen, ann$cis, ann$int, sep="+")
	ann$call <- str_replace_all(ann$call,"NA","")
	ann$call <- str_replace_all(ann$call,"\\+\\+","\\+")
	ann$call <- str_replace_all(ann$call,"^\\+","")
	ann$call <- str_replace_all(ann$call,"\\+$","")

	ann$gen <- NULL
	ann$int <- NULL
	ann$cis <- NULL

	if(gene.names==TRUE)
	{
		print("Pulling gene symbols for overlaps: Genic")
		ann$genic.genes <- getGeneList(input.gr, genic.full.gr, kg)
		print("Pulling gene symbols for overlaps: Flanking")
		ann$flank.genes <- getGeneListForFlanks(input.gr, flank.gr, flank.full.gr, kg)
		print("Pulling gene symbols for overlaps: Nearest")
		ann <- cbind(ann, getNearestGene(input.gr, genic.full.gr, kg))
	}
	ann$url <- getBrowserURLs(input.gr, genome)

	# return the annotated columns
	cbind(query, ann)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Annotate regions in non-mutually exclusive gene model categories
#'
#' Given a set of query regions in a data.frame containing the columns "chr", "start", and "end", returns a new data.frame with the same number of rows as the query and adds annotation columns for the percent of each query region that overlaps with a type of gene model region (5' UTR, 3' UTR, exon, intron, upstream flank, downstream flank). Since there are many overlapping isoforms of genes and instances of overlapping separate genes, the percentages are not mutually exclusive. A "call" column creates discrete categories for each combination of overlaps observed.
#' @param query A data.frame of regions to annotate. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. All additional columns will be retained in the output object.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly.
#' @param flank.bp The size in base pairs for the flanking category.
#' @return A data.frame with all rows and columns of the query data.frame with the annotation columns added.
#' @export
annotateGeneModel <- function(query, genome, cachedir=NULL, flank.bp=1000)
{
	# Check that query has the needed columns
	if(sum(c("chr", "start", "end") %in% colnames(query))!=3)
	{
		stop("Could not find columns named \"chr\", \"start\", and \"end\" in query data.frame")
	}

	print("Generating annotation regions from UCSC knownGene table")

	# Load the UCSC tables we will need
	knownGene <- getUCSCTable("knownGene", genome, cachedir)
	kgXref <- getUCSCTable("kgXref", genome, cachedir)
	knownCanonical <- getUCSCTable("knownCanonical", genome, cachedir)
	chromInfo <- getUCSCTable("chromInfo", genome, cachedir)

	# join in gene symbols
	names(kgXref)[1] <- "name"
	kg <- join(knownGene, kgXref, by="name", type="left", match="first")

	# join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	kg <- join(kg,kgCanon.sub, by="name", type="left", match="first")
	kg[is.na(kg$canonical),]$canonical <- 0

	# add row number field for future rejoins
	kg$srow <- 1:nrow(kg)

	# add row number field to input for future rejoins
	query$qrow <- 1:nrow(query)

	# Turn input data into GRanges
	input.gr <- GRanges(seqnames=query$chr, ranges=IRanges(start=query$start, end=query$end), qrow=query$qrow)
	seqlengths(input.gr) <- chromInfo[match(names(seqlengths(input.gr)), chromInfo$chr),]$size

	# Make genic GRanges
	genic.gr <- with(kg, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1,end=txEnd)))
	genic.gr <- reduce(genic.gr)
	seqlengths(genic.gr) <- chromInfo[match(names(seqlengths(genic.gr)), chromInfo$chr),]$size

	# Make 5' and 3' Flanks
	kg.proms.p <- kg[kg$strand=="+",]
	kg.proms.m <- kg[kg$strand=="-",]

	genic.p.gr <- with(kg.proms.p, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1, end=txEnd)))
	seqlengths(genic.p.gr) <- chromInfo[match(names(seqlengths(genic.p.gr)), chromInfo$chr),]$size
	genic.m.gr <- with(kg.proms.m, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1, end=txEnd)))
	seqlengths(genic.m.gr) <- chromInfo[match(names(seqlengths(genic.m.gr)), chromInfo$chr),]$size

	flank.us.p <- suppressWarnings(flank(genic.p.gr, flank.bp, start=T))
	flank.ds.p <- suppressWarnings(flank(genic.p.gr, flank.bp, start=F))

	flank.us.m <- suppressWarnings(flank(genic.m.gr, flank.bp, start=F))
	flank.ds.m <- suppressWarnings(flank(genic.m.gr, flank.bp, start=T))

	flank.us.gr <- suppressWarnings(reduce(c(flank.us.p, flank.us.m)))
	flank.ds.gr <- suppressWarnings(reduce(c(flank.ds.p, flank.ds.m)))

	# Subtract any overlaps with genic ranges
	flank.us.gr <- setdiff(flank.us.gr, genic.gr)
	flank.ds.gr <- setdiff(flank.ds.gr, genic.gr)

	# Make intergenic range
	intergenic.gr <- GRanges(seqnames=chromInfo$chrom, ranges=IRanges(start=1, end=chromInfo$size))
	seqlengths(intergenic.gr) <- chromInfo[match(names(seqlengths(intergenic.gr)), chromInfo$chr),]$size
	intergenic.gr <- setdiff(intergenic.gr, genic.gr)
	intergenic.gr <- setdiff(intergenic.gr, flank.us.gr)
	intergenic.gr <- setdiff(intergenic.gr, flank.ds.gr)

	# Make exons
	exonstarts <- str_replace(kg$exonStarts,",$","")
	exonstartvec <- unlist(str_split(exonstarts,","))
	exonstartvec <- as.numeric(exonstartvec)+1

	exonends <- str_replace(kg$exonEnds,",$","")
	exonendvec <- as.numeric(unlist(str_split(exonends,",")))

	exonchrs <- rep(kg$chrom, kg$exonCount)

	exon.gr <- GRanges(seqnames=exonchrs, ranges=IRanges(start=exonstartvec,end=exonendvec))

	# Make introns

	## G/IRANGES method
	#mat <- cbind(exonstartvec, exonendvec)

	#tron <- function(x)
	#{
		#gr <- GRanges(seqnames=x[,3], ranges=IRanges(start=x[,1], end=x[,2]))
		#gaps(IRanges(start=x[,1], end=x[,2]))
		#gaps(gr)
	#}
	#tmp <- by(mat, INDICES=rep(paste(kg$name, kg$chrom, sep=":"),kg$exonCount), FUN=tron)
	#tmp2 <- unlist(IRangesList(tmp))
	#tmp3 <- GRanges(seqnames=str_split_fixed(names(tmp2),":",2)[,2], ranges=tmp2)
	##

	## LOOP method
	# We expect to find the number of exons minus 1 for each gene
	print("Computing introns, please wait...")
	intron.mat <- matrix(nrow=sum(kg$exonCount-1), ncol=5)
	r <- 0 #row counter, so we can insert into prealloc matrix for speed boost
	for(j in 1:nrow(kg))
	{
		x <- kg[j,]
		es <- unlist(str_split(x$exonStart,","))
		ee <- unlist(str_split(x$exonEnd,","))

		if(x$exonCount > 1)
		{
			for(i in 1:(x$exonCount-1))
			{
				r <- r+1
				intron.mat[r,1] <- j
				intron.mat[r,2] <- i
				intron.mat[r,3] <- 0
				intron.mat[r,4] <- as.numeric(ee[i])
				intron.mat[r,5] <- as.numeric(es[i+1])

				#new <- c(x$name, paste("Intron",i,sep=" "), x$chrom, ee[i],es[i+1])
				#intron.df <- rbind(intron.df, new)
				#intron.mat[j+1,] <- new
				#print(paste(x$name, paste("Intron",i,sep=" "), x$chrom, ee[i],es[i+1],sep=","))
			}
		}
		if((j %% 10000)==0)
		{
			print(paste("Done with ",j, " genes, found ", r, " introns",sep=""))
		}
	}

	intron.gr <- GRanges(seqnames=kg[intron.mat[,1],]$chrom, ranges=IRanges(start=intron.mat[,4], end=intron.mat[,5]))

	# Make 5' UTR
	# distance between txStart and cdsStart for + genes
	# distance between txEnd and cdsEnd for - genes
	kg.p <- kg[(kg$strand=="+")&(kg$txStart!=kg$cdsStart)&(kg$cdsEnd!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$txEnd!=kg$cdsEnd)&(kg$cdsEnd!=kg$cdsStart),]

	utr5.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chrom, ranges=IRanges(start=txStart+1, end=cdsStart+1))),with(kg.m, GRanges(seqnames=kg.m$chrom, ranges=IRanges(start=cdsEnd, end=txEnd)))))

	# Make 3' UTR
	kg.p <- kg[(kg$strand=="+")&(kg$txEnd!=kg$cdsEnd)&(kg$cdsEnd!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$txStart!=kg$cdsStart)&(kg$cdsEnd!=kg$cdsStart),]

	utr3.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chrom, ranges=IRanges(start=cdsEnd, end=txEnd))),with(kg.m, GRanges(seqnames=kg.m$chrom, ranges=IRanges(start=txStart+1, end=cdsStart+1)))))

	# Make mutually exclusive if option is set (?)

	# write BED files to make sure we have the sets correct
	# writeBEDFromGRanges(genic.gr, file="genic.bed")
	# writeBEDFromGRanges(intergenic.gr, file="intergenic.bed")
	# writeBEDFromGRanges(flank.us.gr, file="flank.us.bed")
	# writeBEDFromGRanges(flank.ds.gr, file="flank.ds.bed")
	# writeBEDFromGRanges(exon.gr, file="exon.bed")
	# writeBEDFromGRanges(intron.gr, file="intron.bed")
	# writeBEDFromGRanges(utr5.gr, file="utr5.bed")
	# writeBEDFromGRanges(utr3.gr, file="utr3.bed")

	#start(reduce(c(genic.gr, intergenic.gr, cis.gr)))
	#intersect(genic.gr, intergenic.gr)
	#intersect(genic.gr, cis.gr)
	#intersect(cis.gr, intergenic.gr)

	# annotate with %
	print("Computing overlaps with gene model range sets")
	ann <- data.frame(intergenic.per=calcPercentOverlap(input.gr, reduce(intergenic.gr)), utr5.per=calcPercentOverlap(input.gr, reduce(utr5.gr)), utr3.per=calcPercentOverlap(input.gr, reduce(utr3.gr)), exon.per=calcPercentOverlap(input.gr, reduce(exon.gr)), intron.per=calcPercentOverlap(input.gr, reduce(intron.gr)), flank.us.per=calcPercentOverlap(input.gr, reduce(flank.us.gr)), flank.ds.per=calcPercentOverlap(input.gr, reduce(flank.ds.gr)))

	# call categories
	print("Calling categories")
	ann$itr <- NA
	ann[ann$intergenic.per>0,]$itr <- "Intergenic"

	ann$utr5 <- NA
	ann[ann$utr5.per>0,]$utr5 <- "5' UTR"

	ann$utr3 <- NA
	ann[ann$utr3.per>0,]$utr3 <- "3' UTR"

	ann$exo <- NA
	ann[ann$exon.per>0,]$exo <- "Exon"

	ann$int <- NA
	ann[ann$intron.per>0,]$int <- "Intron"

	ann$flu <- NA
	ann[ann$flank.us.per>0,]$flu <- "5' Flank"

	ann$fld <- NA
	ann[ann$flank.ds.per>0,]$fld <- "3' Flank"

	ann$call <- paste(ann$itr, ann$utr5, ann$utr3, ann$exo, ann$int, ann$flu, ann$fld, sep="+")
	ann$call <- str_replace_all(ann$call,"NA","")
	ann$call <- str_replace_all(ann$call,"\\++","\\+")
	ann$call <- str_replace_all(ann$call,"^\\+","")
	ann$call <- str_replace_all(ann$call,"\\+$","")

	ann$itr <- NULL
	ann$utr5 <- NULL
	ann$utr3 <- NULL
	ann$exo <- NULL
	ann$int <- NULL
	ann$flu <- NULL
	ann$fld <- NULL

	# add browser links if this was requested
	ann$url <- getBrowserURLs(input.gr, genome)

	# return the annotated columns
	cbind(query, ann)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Report overlap details with one row per pair of overlapping query and gene isoform genomic region
#'
#' Given a set of query regions in a data.frame containing the columns "chr", "start", and "end", returns a new data.frame with one row for each overlapping pair of query region and gene isoform in the UCSC knownGene table. The columns report the percent of the query region that overlaps with gene model regions for that specific isoform (5' UTR, 3' UTR, exons, introns, flanks).
#' @param query A data.frame of regions to annotate. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. All additional columns will be retained in the output object.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly.
#' @param flank.bp The size in base pairs for the flanking category.
#' @return A data.frame with all rows and columns of the query data.frame with the annotation columns added.
#' @export
reportGeneModel <- function(query, genome, cachedir=NULL, flank.bp=1000)
{
	# Check that query has the needed columns
	if(sum(c("chr", "start", "end") %in% colnames(query))!=3)
	{
		stop("Could not find columns named \"chr\", \"start\", and \"end\" in query data.frame")
	}
	chrs <- query$chr
	starts <- query$start
	ends <- query$end

	print("Generating annotation regions from UCSC knownGene table")

	# Load the UCSC tables we will need
	knownGene <- getUCSCTable("knownGene", genome, cachedir)
	kgXref <- getUCSCTable("kgXref", genome, cachedir)
	knownCanonical <- getUCSCTable("knownCanonical", genome, cachedir)
	chromInfo <- getUCSCTable("chromInfo", genome, cachedir)

	# join in gene symbols
	names(kgXref)[1] <- "name"
	kg <- join(knownGene, kgXref, by="name", type="left", match="first")

	# join in canonical annotation (largest coding seq from nearest cluster)
	kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	kg <- join(kg,kgCanon.sub, by="name", type="left", match="first")
	kg[is.na(kg$canonical),]$canonical <- 0

	# add row number field for future rejoins
	kg$srow <- 1:nrow(kg)

	# add row number field to input for future rejoins
	query$qrow <- 1:nrow(query)

	# Turn input data into GRanges
	input.gr <- GRanges(seqnames=query$chr, ranges=IRanges(start=query$start, end=query$end), qrow=query$qrow)
	seqlengths(input.gr) <- chromInfo[match(names(seqlengths(input.gr)), chromInfo$chr),]$size

	# Compute genic overlaps to establish our pairs as a df
	# Then we will join back in each overlap type by matching

	# Make genic GRanges
	kg$srow <- 1:nrow(kg)
	genic.gr <- with(kg, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1,end=txEnd), srow=srow))
	#genic.gr <- reduce(genic.gr)
	seqlengths(genic.gr) <- chromInfo[match(names(seqlengths(genic.gr)), chromInfo$chr),]$size

	fo <- suppressWarnings(findOverlaps(input.gr, genic.gr))
	foq <- queryHits(fo)
	fos <- subjectHits(fo)

	# going to ignore strand for now, so user will have to look at the strand col and know that things like exon numbers are reversed on the "-" strand, else we could also swap after the annotation is done - MAYBE

	qrow=input.gr[foq]$qrow
	srow=genic.gr[fos]$srow
	query.chr=chrs[foq]
	query.start=starts[foq]
	query.end=ends[foq]
	gene.symbol=kg[fos,]$geneSymbol
	isoform.id=kg[fos,]$name
	refseq.id=kg[fos,]$refseq
	isoform.chr=kg[fos,]$chrom
	isoform.start=kg[fos,]$txStart+1
	isoform.end=kg[fos,]$txEnd
	isoform.strand=kg[fos,]$strand
	overlap.bp=suppressWarnings(calcPercentOverlap(input.gr, genic.gr, sum.all=FALSE, report.bp=TRUE))
	query.overlap.per=suppressWarnings(calcPercentOverlap(input.gr, genic.gr, sum.all=FALSE))
	isoform.overlap.per=suppressWarnings(calcPercentOverlap(genic.gr, input.gr, sum.all=FALSE))

	ann <- data.frame(qrow=qrow, srow=srow, query.chr=query.chr, query.start=query.start, query.end=query.end, gene.symbol=gene.symbol, isoform.id=isoform.id, refseq.id=refseq.id, isoform.chr=isoform.chr, isoform.start=isoform.start, isoform.end=isoform.end, isoform.strand=isoform.strand, overlap.bp=overlap.bp, query.overlap.per=query.overlap.per, isoform.overlap.per=isoform.overlap.per)
	#gene.desc=kg[fos,]$description
	print("Initial Intersect Complete")
	# Make 5' and 3' Flanks
	kg.proms.p <- kg[kg$strand=="+",]
	kg.proms.m <- kg[kg$strand=="-",]

	# Going to carry the srow from kg to map matches back to their source rows
	genic.p.gr <- with(kg.proms.p, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1, end=txEnd),srow=srow))
	seqlengths(genic.p.gr) <- chromInfo[match(names(seqlengths(genic.p.gr)), chromInfo$chr),]$size
	genic.m.gr <- with(kg.proms.m, GRanges(seqnames=chrom, ranges=IRanges(start=txStart+1, end=txEnd),srow=srow))
	seqlengths(genic.m.gr) <- chromInfo[match(names(seqlengths(genic.m.gr)), chromInfo$chr),]$size

	flank.us.p <- suppressWarnings(flank(genic.p.gr, flank.bp, start=T))
	flank.ds.p <- suppressWarnings(flank(genic.p.gr, flank.bp, start=F))

	flank.us.m <- suppressWarnings(flank(genic.m.gr, flank.bp, start=F))
	flank.ds.m <- suppressWarnings(flank(genic.m.gr, flank.bp, start=T))

	flank.us.gr <- suppressWarnings(c(flank.us.p, flank.us.m))
	flank.ds.gr <- suppressWarnings(c(flank.ds.p, flank.ds.m))

	# Subtract any overlaps with genic ranges
	# NOT going to subtract this time - flanks MAY include overlaps with other genes on this format (it is too complicated to deal with subtracting and having flanks potentially cut up and having to map them back. I need to know, one gene, one flank for the future intersect code.)
	#flank.us.gr2 <- setdiff(flank.us.gr, genic.gr)
	#flank.ds.gr2 <- setdiff(flank.ds.gr, genic.gr)

	flank.us.o <- calcOverlapForReport(input.gr, flank.us.gr, all=FALSE)
	colnames(flank.us.o) <- c("qrow", "srow", "5' Flank")
	ann <- join(ann, flank.us.o, type="left", by=c("qrow", "srow"))
	ann[is.na(ann$"5' Flank"),]$"5' Flank" <- 0


	# Make 5' UTR
	# distance between txStart and cdsStart for + genes
	# distance between txEnd and cdsEnd for - genes
	kg.p <- kg[(kg$strand=="+")&(kg$txStart!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$txEnd!=kg$cdsEnd),]

	utr5.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chrom, ranges=IRanges(start=txStart+1, end=cdsStart+1), srow=srow)),with(kg.m, GRanges(seqnames=kg.m$chrom, ranges=IRanges(start=cdsEnd, end=txEnd), srow=srow))))

	utr5.o <- calcOverlapForReport(input.gr, utr5.gr, all=FALSE)
	colnames(utr5.o) <- c("qrow", "srow", "5' UTR")
	ann <- join(ann, utr5.o, type="left", by=c("qrow", "srow"))
	ann[is.na(ann$"5' UTR"),]$"5' UTR" <- 0




	print("Flanks and UTRs Computed")



	# Make intergenic range
	#intergenic.gr <- GRanges(seqnames=chromInfo$chrom, ranges=IRanges(start=1, end=chromInfo$size))
	#seqlengths(intergenic.gr) <- chromInfo[match(names(seqlengths(intergenic.gr)), chromInfo$chr),]$size
	#intergenic.gr <- setdiff(intergenic.gr, genic.gr)
	#intergenic.gr <- setdiff(intergenic.gr, flank.us.gr)
	#intergenic.gr <- setdiff(intergenic.gr, flank.ds.gr)

	# Make exons
	exonstarts <- str_replace(kg$exonStarts,",$","")
	exonstartvec <- unlist(str_split(exonstarts,","))
	exonstartvec <- as.numeric(exonstartvec)+1

	exonends <- str_replace(kg$exonEnds,",$","")
	exonendvec <- as.numeric(unlist(str_split(exonends,",")))

	exonchrs <- rep(kg$chrom, kg$exonCount)

	exonsrows <- rep(kg$srow, kg$exonCount)

	exonnums <- unlist(sapply(kg$exonCount, FUN=function(x) seq(1, x)))

	exon.gr <- GRanges(seqnames=exonchrs, ranges=IRanges(start=exonstartvec,end=exonendvec), srow=exonsrows, num=exonnums)

	# similar to the single-time regions, we just need to keep track of the intron and exon numbers and carry these through so we can position our list properly (and remember to flip it around for "-" genes!)
	exon.o <- calcOverlapForReportExons(input.gr, exon.gr)

	#order(exon.o$qrow, exon.o$srow, exon.o$num)

	#procex <- function(x)
	#{
		#E1(25%), I1(25%), I2(0%), E2(0%)
	#	print(paste("E", paste(x$num, x$per, collapse=" "), sep="", collapse=","))
	#}
	#down <- by(exon.o, INDICES=paste(exon.o$qrow, exon.o$srow,sep="-"), FUN=procex)

	# Make introns

	## G/IRANGES method
	#mat <- cbind(exonstartvec, exonendvec)

	#tron <- function(x)
	#{
		#gr <- GRanges(seqnames=x[,3], ranges=IRanges(start=x[,1], end=x[,2]))
		#gaps(IRanges(start=x[,1], end=x[,2]))
		#gaps(gr)
	#}
	#tmp <- by(mat, INDICES=rep(paste(kg$name, kg$chrom, sep=":"),kg$exonCount), FUN=tron)
	#tmp2 <- unlist(IRangesList(tmp))
	#tmp3 <- GRanges(seqnames=str_split_fixed(names(tmp2),":",2)[,2], ranges=tmp2)
	##

	## LOOP method
	# We expect to find the number of exons minus 1 for each gene
	print("Computing introns, please wait...")
	intron.mat <- matrix(nrow=sum(kg$exonCount-1), ncol=5)
	r <- 0 #row counter, so we can insert into prealloc matrix for speed boost
	for(j in 1:nrow(kg))
	{
		x <- kg[j,]
		es <- unlist(str_split(x$exonStart,","))
		ee <- unlist(str_split(x$exonEnd,","))

		if(x$exonCount > 1)
		{
			for(i in 1:(x$exonCount-1))
			{
				r <- r+1
				intron.mat[r,1] <- j
				intron.mat[r,2] <- i
				intron.mat[r,3] <- 0
				intron.mat[r,4] <- as.numeric(ee[i])
				intron.mat[r,5] <- as.numeric(es[i+1])

				#new <- c(x$name, paste("Intron",i,sep=" "), x$chrom, ee[i],es[i+1])
				#intron.df <- rbind(intron.df, new)
				#intron.mat[j+1,] <- new
				#print(paste(x$name, paste("Intron",i,sep=" "), x$chrom, ee[i],es[i+1],sep=","))
			}
		}
		if((j %% 10000)==0)
		{
			print(paste("Done with ",j, " genes, found ", r, " introns",sep=""))
		}
	}

	intron.gr <- GRanges(seqnames=kg[intron.mat[,1],]$chrom, ranges=IRanges(start=intron.mat[,4]+1, end=intron.mat[,5]), srow=kg[intron.mat[,1],]$srow, num=intron.mat[,2])

	intron.o <- calcOverlapForReportExons(input.gr, intron.gr)

	exon.o$type <- "E"
	intron.o$type <- "I"

	all.o <- rbind(exon.o, intron.o)

	# Make a DF with all possible introns and exon numbers
	print("Generating exon/intron overlaps")
	kgxs <- rbind(data.frame(srow=exon.gr$srow, num=exon.gr$num, type="E"), data.frame(srow=intron.gr$srow, num=intron.gr$num, type="I"))

	pairs <- data.frame(qrow=ann$qrow, srow=ann$srow)

	pairs.dt <- data.table(pairs)
	setkey(pairs.dt, qrow, srow)
	kgxs.dt <- data.table(kgxs)
	setkey(kgxs.dt, srow)

	kgx.dt <- merge(pairs.dt, kgxs.dt, by="srow", allow.cartesian=TRUE)

	# Join our percents in
	all.o.dt <- data.table(all.o)
	setkey(all.o.dt, qrow, srow, num, type)
	ts.dt <- merge(kgx.dt, all.o.dt, by=c("qrow","srow","num","type"), all.x=TRUE, all.y=FALSE)

	#ts <- join(kgx, all.o, by=c("qrow","srow","type"))

	# Assign the rest as 0%
	if(sum(is.na(ts.dt$per))>0)
	{
		ts.dt[is.na(ts.dt$per),]$per <- 0
	}

	print("Processing exon/intron overlap diagrams")

	# Now use this to make the format strings
	i<-0
	procex <- function(x)
	{
		i<<-i+1
		if((i %% 10000)==0)
		{
			print(paste("Made exon/intron strings for ",i, " overlaps",sep=""))
		}
		#E1(25%), I1(25%), I2(0%), E2(0%)
		ourstrand <- kg[x[1,]$srow,]$strand
		if(ourstrand=="-")
		{
			# Flip the Order if we are "-"
			x[x$type=="E",]$num <- rev(x[x$type=="E",]$num)
			x[x$type=="I",]$num <- rev(x[x$type=="I",]$num)
			x <- x[order(x$num, x$type),]
		}

		x2 <- x[order(x$num),]
		x2$string <- paste(x2$type, x2$num, " (", x2$per,")", sep="")
		paste(x2$string, collapse=", ")
	}
	ts.dt$pair <- paste(ts.dt$qrow, ts.dt$srow,sep="-")
	down2 <- ts.dt[,procex(data.frame(qrow, srow, num, type, per)),by=pair]

	##NEXT
	down <- data.frame(down2)
	ann$pair <- paste(ann$qrow, ann$srow,sep="-")
	ann <- join(ann, down2, by="pair", match="first")
	ann$pair <- NULL
	ann$"Introns and Exons" <- ann$V1
	ann$V1 <- NULL

	#ann$"Exons and Introns" <- down
	# want a function where given a subset of this DF, return string formatted as:
	# E1 (X%), I1 (X%), E2 (X%), E3 (X%)
	# But also want the zeros to show up


	# Make 3' UTR
	kg.p <- kg[(kg$strand=="+")&(kg$txEnd!=kg$cdsEnd)&(kg$cdsEnd!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$txStart!=kg$cdsStart)&(kg$cdsEnd!=kg$cdsStart),]

	utr3.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chrom, ranges=IRanges(start=cdsEnd, end=txEnd), srow=srow)),with(kg.m, GRanges(seqnames=kg.m$chrom, ranges=IRanges(start=txStart+1, end=cdsStart+1), srow=srow))))



	utr3.o <- calcOverlapForReport(input.gr, utr3.gr, all=FALSE)
	colnames(utr3.o) <- c("qrow", "srow", "3' UTR")
	ann <- join(ann, utr3.o, type="left", by=c("qrow", "srow"))
	ann[is.na(ann$"3' UTR"),]$"3' UTR" <- 0

	flank.ds.o <- calcOverlapForReport(input.gr, flank.ds.gr, all=FALSE)
	colnames(flank.ds.o) <- c("qrow", "srow", "3' Flank")
	ann <- join(ann, flank.ds.o, type="left", by=c("qrow", "srow"))
	ann[is.na(ann$"3' Flank"),]$"3' Flank" <- 0


	# ## PATHWAY STUFF
	# bioCycPathway <- getUCSCTable("bioCycPathway", genome, cachedir)
	# bioCycMapDesc <- getUCSCTable("bioCycMapDesc", genome, cachedir)

	# #biocyc
	# ann$kgID <- ann$isoform.id
	# ann <- join(ann, bioCycPathway, by="kgID", match="first")
	# #unique(y$mapID)
	# ann <- join(ann, bioCycMapDesc, by="mapID", match="first")
	# ann$kgID <- NULL
	# ann$geneID <- NULL
	# ann$mapID <- NULL
	# ann$bioCycMapDesc <- ann$description
	# ann[is.na(ann$bioCycMapDesc),]$bioCycMapDesc <- ""
	# ann$description <- NULL 

	# #kegg
	# keggPathway <- getUCSCTable("keggPathway", genome, cachedir)
	# keggMapDesc <- getUCSCTable("keggMapDesc", genome, cachedir)
	# ann$kgID <- ann$isoform.id
	# ann <- join(ann, keggPathway, by="kgID", match="first")
	# #unique(y$mapID)
	# ann <- join(ann, keggMapDesc, by="mapID", match="first")
	# ann$kgID <- NULL
	# ann$locusID <- NULL
	# ann$mapID <- NULL
	# ann$keggMapDesc <- ann$description
	# ann[is.na(ann$keggMapDesc),]$keggMapDesc <- ""
	# ann$description <- NULL 

	#go
	#ann$proteinID <- kg[fos,]$proteinID

	#goaPart <- getUCSCTable("goaPart", "go", cachedir)
	#term <- getUCSCTable("term", "go", cachedir)

	#goaPart$proteinID <- goaPart$dbObjectId
	#ann <- join(ann, goaPart, by="proteinID", match="first")

	#gene desc
	#ann$gene.desc <- kg[fos,]$description
	#write.csv(ann, "report.try4.csv", row.names=FALSE)




	# Make mutually exclusive if option is set (?)

	# write BED files to make sure we have the sets correct
	# writeBEDFromGRanges(genic.gr, file="genic.bed")
	# writeBEDFromGRanges(intergenic.gr, file="intergenic.bed")
	# writeBEDFromGRanges(flank.us.gr, file="flank.us.bed")
	# writeBEDFromGRanges(flank.ds.gr, file="flank.ds.bed")
	# writeBEDFromGRanges(exon.gr, file="exon.bed")
	# writeBEDFromGRanges(intron.gr, file="intron.bed")
	# writeBEDFromGRanges(utr5.gr, file="utr5.bed")
	# writeBEDFromGRanges(utr3.gr, file="utr3.bed")

	#start(reduce(c(genic.gr, intergenic.gr, cis.gr)))
	#intersect(genic.gr, intergenic.gr)
	#intersect(genic.gr, cis.gr)
	#intersect(cis.gr, intergenic.gr)

	# annotate with %
	#ann <- data.frame(intergenic.per=calcPercentOverlap(input.gr, intergenic.gr), utr5.per=calcPercentOverlap(input.gr, utr5.gr), utr3.per=calcPercentOverlap(input.gr, utr3.gr), exon.per=calcPercentOverlap(input.gr, exon.gr), intron.per=calcPercentOverlap(input.gr, intron.gr), flank.us.per=calcPercentOverlap(input.gr, flank.us.gr), flank.ds.per=calcPercentOverlap(input.gr, flank.ds.gr))

	# call categories
	# ann$itr <- NA
	# ann[ann$intergenic.per>0,]$itr <- "Intergenic"

	# ann$utr5 <- NA
	# ann[ann$utr5.per>0,]$utr5 <- "5' UTR"

	# ann$utr3 <- NA
	# ann[ann$utr3.per>0,]$utr3 <- "3' UTR"

	# ann$exo <- NA
	# ann[ann$exon.per>0,]$exo <- "Exon"

	# ann$int <- NA
	# ann[ann$intron.per>0,]$int <- "Intron"

	# ann$flu <- NA
	# ann[ann$flank.us.per>0,]$flu <- "5' Flank"

	# ann$fld <- NA
	# ann[ann$flank.ds.per>0,]$fld <- "3' Flank"

	# ann$call <- paste(ann$itr, ann$utr5, ann$utr3, ann$exo, ann$int, ann$flu, ann$fld, sep="+")
	# ann$call <- str_replace_all(ann$call,"NA","")
	# ann$call <- str_replace_all(ann$call,"\\++","\\+")
	# ann$call <- str_replace_all(ann$call,"^\\+","")
	# ann$call <- str_replace_all(ann$call,"\\+$","")

	# ann$itr <- NULL
	# ann$utr5 <- NULL
	# ann$utr3 <- NULL
	# ann$exo <- NULL
	# ann$int <- NULL
	# ann$flu <- NULL
	# ann$fld <- NULL

	# add noncoding flag
	ann$noncoding <- kg[fos,]$cdsStart==kg[fos,]$cdsEnd

	# add browser links if this was requested
	ann$url <- getBrowserURLs(GRanges(seqnames=ann$query.chr, ranges=IRanges(start=ann$query.start,end=ann$query.end)), genome)

	# return the annotated columns
	cbind(query[foq,], ann)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Annotate percent of query that overlaps with any regions in features, or report unique IDs from features for a specified column that overlap with query
#'
#' Given a set of query regions in a data.frame containing the columns "chr", "start", and "end", returns a vector in the same length and order as the rows in query with either the percent of each query that overlap with any regions in features (default), or, if getid is given as a string, the function will report a comma separated string of IDs that are involved. Raw UCSC tables can be passed directly if ucsc=TRUE is given, and this option will adjust for 0-based start coordinates and UCSC's column names for the coordinates.
#' @param query A data.frame of regions to annotate. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. Columns will not be retained in output, but the output vector will match query row for row.
#' @param features A data.frame of regions to overlap query against. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. Or, must be a raw UCSC table from getUCSCTable() if ucsc=TRUE is also given.
#' @param getid If NULL, return percentage overlaps. If a string, the name of a column in features to pull overlapping IDs or other data fields from.
#' @param ucsc Set to TRUE if the table is a raw UCSC table fetched with getUCSCTable()
#' @return A vector of the requested overlap data, matching row for row with query.
#' @export
annotateFeatures <- function(query, features, getid=NULL, ucsc=FALSE)
{
	# Check that query has the needed columns
	if(sum(c("chr", "start", "end") %in% colnames(query))!=3)
	{
		stop("Could not find columns named \"chr\", \"start\", and \"end\" in query data.frame")
	}
	if(ucsc==TRUE)
	{
		if(sum(c("chrom", "chromStart", "chromEnd") %in% colnames(features))==3)
		{
			print("Treating as UCSC raw table. Converting start coordinates from 0-based to 1-based.")
			features.chr <- features$chrom
			features.start <- features$chromStart + 1
			features.end <- features$chromEnd
		} else
		{
			stop("Could not find UCSC columns chrom, chromStart, and chromEnd in features dataframe.")
		}
	} else
	{
		if(sum(c("chr", "start", "end") %in% colnames(features))==3)
		{
			features.chr <- features$chr
			features.start <- features$start
			features.end <- features$end
		} else
		{
			stop("Could not find chr, start, and end columns in features dataframe.")
		}
	}

	if((!is.null(getid)))
	{
		if((!(getid %in% colnames(features))))
		{
			stop("Could not find get.id column in features dataframe.")
		}
	}

	# Make input GRanges
	input.gr <- GRanges(seqnames=query$chr, ranges=IRanges(start=query$start, end=query$end), qrow=1:nrow(query))

	# Make feature GRanges
	features.gr <- GRanges(seqnames=features.chr, ranges=IRanges(start=features.start, end=features.end), srow=1:nrow(features))
	features.full.gr <- features.gr
	features.gr <- reduce(features.gr)

	getOverlapIds <- function(query.gr, subject.gr, features, getid)
	{
		# Perform overlaps
		overlaps.fo <- findOverlaps(query.gr, subject.gr)

		# Make DF with gene symbol for each overlap
		genes.df <- data.frame(qrow=queryHits(overlaps.fo), srow=subjectHits(overlaps.fo))
		genes.df$gene <- features[genes.df$srow,colnames(features) %in% getid]

		# Condense into comma-separated lists, removing duplicate symbols
		genes.dt <- data.table(genes.df)
		setkey(genes.dt, qrow, srow)
		i<-0
		commaGenes <- function(x)
		{
			i<<-i+1
			if((i %% 10000)==0)
			{
				print(paste("Done with feature names for ",i, " regions",sep=""))
			}
			genes <- unique(x)
			genes <- paste(genes, collapse=", ")
		}
		g <- genes.dt[,commaGenes(gene),by=qrow]

		# Add blanks for intergenics
		out <- data.table(qrow=query.gr$qrow)
		out <- merge(out, g, by="qrow", all.x=TRUE, all.y=FALSE)
		setkey(out, qrow)
		if(sum(is.na(out$V1)) > 0)
		{
			out[is.na(V1),]$V1 <- ""
		}
		out$V1
	}

	# Report overlap
	if(is.null(getid))
	{
		calcPercentOverlap(input.gr, features.gr, sum.all=TRUE, report.bp=FALSE)
	} else
	{
		getOverlapIds(input.gr, features.full.gr, features, getid)
	} 
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#' Report one row per pair of overlapping query and feature range. Retains all columns from both query and features.
#'
#' Given a set of query regions in a data.frame containing the columns "chr", "start", and "end", returns a data.frame with all columns from both query and features for each pair of overlapping ranges from each group. Raw UCSC tables can be passed directly if ucsc=TRUE is given, and this option will adjust for 0-based start coordinates and UCSC's column names for the coordinates.
#' @param query A data.frame of regions to annotate. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. Columns will not be retained in output, but the output vector will match query row for row.
#' @param features A data.frame of regions to overlap query against. Must contain the columns "chr", "start", "end", and the "start" coordinates must be 1-based. Or, must be a raw UCSC table from getUCSCTable() if ucsc=TRUE is also given.
#' @param ucsc Set to TRUE if the table is a raw UCSC table fetched with getUCSCTable()
#' @return A data.frame, where each row is an overlapping pair of query and features regions, and contains data about the overlap and all columns from both data sets.
#' @export
reportFeatures <- function(query, features, ucsc=FALSE)
{
	if(sum(c("chr", "start", "end") %in% colnames(query))!=3)
	{
		stop("Could not find columns named \"chr\", \"start\", and \"end\" in query data.frame")
	}
	if(ucsc==TRUE)
	{
		if(sum(c("chrom", "chromStart", "chromEnd") %in% colnames(features))==3)
		{
			print("Treating as UCSC raw table. Converting start coordinates from 0-based to 1-based.")
			features.chr <- features$chrom
			features.start <- features$chromStart + 1
			features.end <- features$chromEnd

			features$bin <- NULL
			features$chrom <- NULL
			features$chromStart <- NULL
			features$chromEnd <- NULL
		} else
		{
			stop("Could not find UCSC columns chrom, chromStart, and chromEnd in features dataframe.")
		}
	} else
	{
		if(sum(c("chr", "start", "end") %in% colnames(features))==3)
		{
			features.chr <- features$chr
			features.start <- features$start
			features.end <- features$end
		} else
		{
			stop("Could not find chr, start, and end columns in features dataframe.")
		}
	}

	# Make input GRanges
	chrs <- query$chr
	starts <- query$start
	ends <- query$end
	input.gr <- GRanges(seqnames=query$chr, ranges=IRanges(start=query$start, end=query$end), qrow=1:nrow(query))

	# Make feature GRanges
	features.gr <- GRanges(seqnames=features.chr, ranges=IRanges(start=features.start, end=features.end), srow=1:nrow(features))

	overs <- calcOverlapForReport(input.gr, features.gr)

	ann <- data.frame(query.chr=chrs[overs$qrow], query.start=starts[overs$qrow], query.end=ends[overs$qrow], feature.chr=features.chr[overs$srow], feature.start=features.start[overs$srow], feature.end=features.end[overs$srow], overlap.query.per=overs$per, overlap.feature.per=overs$sper, overlap.bp=overs$bp)
	ann <- cbind(query[overs$qrow,], ann, features[overs$srow,])
}
# -----------------------------------------------------------------------------




# =============================================================================

# =============================================================================
# Internal Functions

# -----------------------------------------------------------------------------
# Internal function for ggplot2 themes
ggnice <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(color="black"))
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal function to calculate % overlaps for report format
# special version to retain a srow and return this back to us for mapping back to the query vs genes intersect set
# reports percent of GENE MODEL FEATURE overlapped, not percent of query that overlaps with it
# return DF of srow from original kg and qrow from original input ranges. We'll then use these to join these features back into the original query<->genes overlap table we're already started contructing. The NAs can then simply be cast over to 0.
calcOverlapForReport <- function(query.gr, subject.gr, all=TRUE)
{
	#subject.gr <- reduce(subject.gr)
	overlaps.fo <- findOverlaps(query.gr, subject.gr)

	# compute width of each overlap row
	ranges.fo <- ranges(overlaps.fo, ranges(query.gr), ranges(subject.gr))

	# make df with each query and its overlapping width
	ofoq <- queryHits(overlaps.fo)
	ofos <- subjectHits(overlaps.fo)
	width.df <- data.frame(query=ofoq, subject=ofos, qrow=query.gr[ofoq]$qrow, srow=subject.gr[ofos]$srow, width=width(ranges.fo))
	width.df$per <- round(width.df$width/width(query.gr[width.df$query])*100,digits=2)
	width.df$sper <- round(width.df$width/width(subject.gr[width.df$subject])*100,digits=2)

	# return DF
	if(all==TRUE)
	{
		data.frame(qrow=width.df$qrow, srow=width.df$srow, per=width.df$per, bp=width.df$width, sper=width.df$sper)
	} else
	{
		data.frame(qrow=width.df$qrow, srow=width.df$srow, per=width.df$per)
	}
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal format to report % overlaps for exons/introns (retains exon # in returned DF)
calcOverlapForReportExons <- function(query.gr, subject.gr, sum.all=TRUE, report.bp=FALSE)
{
	#subject.gr <- reduce(subject.gr)
	overlaps.fo <- findOverlaps(query.gr, subject.gr)

	# compute width of each overlap row
	ranges.fo <- ranges(overlaps.fo, ranges(query.gr), ranges(subject.gr))

	# make df with each query and its overlapping width
	ofoq <- queryHits(overlaps.fo)
	ofos <- subjectHits(overlaps.fo)
	width.df <- data.frame(query=ofoq, subject=ofos, qrow=query.gr[ofoq]$qrow, srow=subject.gr[ofos]$srow, num=subject.gr[ofos]$num, width=width(ranges.fo))
	width.df$per <- round(width.df$width/width(query.gr[width.df$query])*100,digits=2)

	# return DF
	data.frame(qrow=width.df$qrow, srow=width.df$srow, num=width.df$num, per=width.df$per)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# internal function to take two GRanges objects and output the % overlap of the query for each of the queries
# expands the output to include 0%s - should match row number of query.gr
calcPercentOverlap <- function(query.gr, subject.gr, sum.all=TRUE, report.bp=FALSE)
{
	#subject.gr <- reduce(subject.gr)
	overlaps.fo <- findOverlaps(query.gr, subject.gr)

	# compute width of each overlap row
	ranges.fo <- ranges(overlaps.fo, ranges(query.gr), ranges(subject.gr))

	# make df with each query and its overlapping width
	width.df <- data.frame(query=queryHits(overlaps.fo), width=width(ranges.fo))

	if(sum.all==TRUE)
	{
		# aggregate down to get the total width for each unique query id
		width.ag <- tapply(width.df$width, width.df$query, FUN=sum)
		width.ag <- data.frame(query=as.numeric(names(width.ag)), width=width.ag)

		# tapply was way faster
		#width.ag <- ddply(width.df, .(query), summarize, width=sum(width))

		# add a column for the total size of the entire query
		width.ag$querySize <- width(query.gr)[width.ag$query]

		# compute what percent of this query is the overlap
		width.ag$percent <- round((width.ag$width / width.ag$querySize)*100, digits=2)

		# join back into a vector for all
		all <- data.frame(query=as.numeric(seq(1, length(query.gr))))

		all.j <- merge(all, width.ag, all.x=TRUE)

		#all.j <- join(all, width.ag, by="query")
		pervec <- all.j$percent
		pervec[is.na(pervec)] <- 0

		# return vector
		pervec
	} else if (report.bp==TRUE)
	{
			pervec <- width.df$width
	} else
	{
			pervec <- round(width.df$width/width(query.gr[width.df$query])*100,digits=2)
	}
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Internal function to dump a GRanges to a BED file (used to double check our ranges visually in UCSC)
writeBEDFromGRanges <- function(gr, file)
{
	fileConn<-file(file)
	writeLines(c(paste("track name=\"",file,"\"",sep="")), fileConn)
	close(fileConn)
	df <- data.frame(chr=as.character(seqnames(gr)),start=start(gr), end=end(gr))
	df <- df[df$chr %in% c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"),]
	write.table(df, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# internal fucntion to take vectors of positions and a genome and output direct URLs to UCSC genome browser
getBrowserURLs <- function(input.gr, genome)
{
	# Give links to UCSC at the position of this region
	#http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr22%3A1-500
	paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=", genome, "&position=", seqnames(input.gr), "%3A", start(input.gr), "-", end(input.gr),sep="")
}
# -----------------------------------------------------------------------------

# =============================================================================