# Basic useful functions

# ====================================================================
# Exported Functions

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
			if(sum(colnames(obj) %in% skipcols)!=ncol(obj))
			{
				values(ret) <- obj[,!(colnames(obj) %in% skipcols),with=F]
			}
		} else
		{
			if(sum(colnames(obj) %in% skipcols)!=ncol(obj))
			{
				values(ret) <- obj[,!(colnames(obj) %in% skipcols)]
			}
		}
		ret
	} else
	{
		# Some bad input, throw error
		stop("GRanges, data.frame, or data.table object required as input")
	}
}
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
#' Write a BED format file from a GenomicRanges object
#'
#' Creates BED file suitable for upload as a custom track to the UCSC genome browser. Note that start coordinates are 0-based in the BED format.
#' @param gr A GenomicRanges object.
#' @param file Filename of the BED file to write.
#' @param name Column name to use for the name field in the BED file (optional)
#' @export
writeBEDFromGRanges <- function(gr, file, name=NULL)
{
	#if(file.exists(file)){file.remove(file)}
	#fileConn<-file(file)
	#writeLines(c(paste("track name=\"",file,"\"",sep="")), fileConn)
	#close(fileConn)
	df <- data.frame(chr=as.character(seqnames(gr)),start=start(gr)-1, end=end(gr))
	if(!is.null(name))
	{
		df$name <- values(gr)[,name]
	}
	#df <- df[df$chr %in% c(sapply(seq(1,22),function(x) paste("chr",x,sep="")),"chrX","chrY"),]
	write.table(df, file=file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=FALSE)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Clean ggplot2 theme
#
# Remove gridlines from ggplot2.
ggnice <- function()
{
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(color="black"))
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Make a data.table from a GRanges or a data.frame
#'
#' Given a data.frame or GRanges, a data.table object will be created. If the input is already a data.table, it is simply returned.
#' @param obj A data.frame or GRanges
#' @return A data.table made from the data in obj.
#' @export
makeDT <- function(obj)
{
	myclass <- class(obj)[1]
	if(myclass=="GRanges")
	{
		# Convert GRanges to DT
		obj.dt <- data.table(as(obj,"data.frame"))
		setnames(obj.dt,"seqnames","chr")
		return(obj.dt)
	} else if(myclass=="data.frame")
	{
		# Convert DF to DT
		return(data.table(obj))
	} else if(myclass=="data.table")
	{
		# Already a DT, just return it back
		return(data.table(obj))
	} else
	{
		# Some bad input, throw error
		stop("GRanges, data.frame, or data.table object required as input")
	}
}
# --------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Add columns to query with distance to nearest subject and subject id(s)
#'
#' @param query Genomic regions to find nearest genes for as a GRanges, data.frame, or data.table with the columns "chr", "start", and "end"
#' @param query Genomic regions to find nearest genes for as a GRanges, data.frame, or data.table with the columns "chr", "start", and "end"
#' @param id Column name of the id field in subject to report as the nearest id(s). In case of ties, a comma separated list will be returned.
#' @param prefix Append this string to names of the added columns
#' @export
addNearest <- function(query,subject,id="name",prefix="subject")
{
	#subject.dt <- as(as(subject,"data.frame"),"data.table")
	query.dt <- makeDT(query)
	subject.dt <- makeDT(subject)

	if(!(id %in% colnames(subject.dt))){stop(paste0("id field \"",id,"\" is not a column in subject"))}

	query.gr <- makeGRanges(query)
	subject.gr <- makeGRanges(subject)

	# Get distances
	dist <- as.data.frame(distanceToNearest(query.gr,subject.gr))$distance

	# Get names
	fo <- data.table(as.data.frame(nearest(query.gr,subject.gr,select="all")))
	fo$name <- subject.dt[fo$subjectHits,][[id]]
	fo <- fo[,list(name2=toString(unique(name))),by=queryHits]
	matched <- data.frame(id=1:length(query.gr),name=NA)
	matched[fo$queryHits,]$name <- fo$name2

	query.dt[,eval(paste(prefix,"nearest",sep="_")):=matched$name]
	query.dt[,eval(paste(prefix,"dist",sep="_")):=dist]
	return(query.dt)
}
# --------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Sort a data.frame, data.table, or GRanges by chr (accounting for mixed string and numeric names), start, end and return a data.table
#'
#' @param obj A data.frame, data.table, or GRanges
#' @export
sortDT <- function(obj)
{
	obj <- makeDT(obj)
	if(sum(c("chr", "start", "end") %in% colnames(obj))!=3){stop("Could not find columns named \"chr\", \"start\", and \"end\" in input data.frame")}
	chrorder <- mixedsort(unique(obj$chr))
	obj$chr <- factor(obj$chr,levels=chrorder)
	obj <- obj[order(obj$chr,obj$start,obj$end),]
	return(obj)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Sort a data.frame, data.table, or GRanges by chr (accounting for mixed string and numeric names), start, end and return a GRanges
#'
#' @param obj A data.frame, data.table, or GRanges
#' @export
sortGRanges <- function(obj)
{
	obj <- makeGRanges(obj)
	chrorder <- mixedsort(unique(seqnames(obj)))
	seqlevels(obj) <- as.character(chrorder)
	obj <- obj[order(start(obj),end(obj))]
	obj <- obj[order(seqnames(obj))]
	return(obj)
}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Add columns with distance to nearest gene and gene symbol(s)
#'
#' @param query Genomic regions to find nearest genes for as a GRanges, data.frame, or data.table with the columns "chr", "start", and "end"
#' @param geneset Select one of "ucsc" for the UCSC Genes (from the knownGene table), "refseq" for RefSeq genes (from the refFlat table), or "ensembl" for the Ensembl genes (from the ensGene table)
#' @param genome UCSC genome name to use (e.g. hg19, mm10)
#' @param cachedir Path where cached UCSC tables are stores
#' @param sync If TRUE, then check if newer versions of UCSC tables are available and download them if so. If FALSE, skip this check. Can be used to freeze data versions in an analysis-specific cachedir for reproducibility.
#' @export
addGenes <- function(query,geneset,genome,cachedir,sync=TRUE)
{
	genes <- getGenes(geneset,genome,cachedir,sync=sync)
	addNearest(query,genes,id="name",prefix=geneset)
}
# --------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' Load table of gene ranges via UCSC Genome Browser tables
#'
#' @param geneset Select one of "ucsc" for the UCSC Genes (from the knownGene table), "refseq" for RefSeq genes (from the refFlat table), or "ensembl" for the Ensembl genes (from the ensGene table)
#' @param genome UCSC genome name to use (e.g. hg19, mm10)
#' @param cachedir Path where cached UCSC tables are stores
#' @param sync If TRUE, then check if newer versions of UCSC tables are available and download them if so. If FALSE, skip this check. Can be used to freeze data versions in an analysis-specific cachedir for reproducibility.
#' @export
getGenes <- function(geneset="ucsc",genome,cachedir=NULL,sync=TRUE)
{
	# Validate geneset
	if(!(geneset %in% c("ucsc","refseq","ensembl"))){stop("geneset must be one of \"ucsc\", \"refseq\", or \"ensembl\"")}

	if(geneset=="ucsc")
	{
		kg <- getUCSCTable("knownGene",genome,cachedir,sync=sync)
		kgx <- suppressWarnings(getUCSCTable("kgXref",genome,cachedir,sync=sync))
		ki <- getUCSCTable("knownIsoforms",genome,cachedir,sync=sync)
		setnames(kg,"name","kgID")
		setnames(ki,"transcript","kgID")
		setkey(kg,kgID)
		setkey(kgx,kgID)
		setkey(ki,kgID)
		kg <- ki[kg,]
		genes <- kgx[kg,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=geneSymbol,gene.id=clusterId,isoform.id=NA,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds,kgID=kgID)]
		genes[,isoform.id:=kgID]
		genes[,kgID:=NULL]
		return(genes)
	} else if(geneset=="refseq")
	{
		rg <- getUCSCTable("refFlat",genome,cachedir,sync=sync)
		genes <- rg[,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=geneName,gene.id=geneName,isoform.id=name,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds)]
		return(genes)
	} else if(geneset=="ensembl")
	{
		eg <- getUCSCTable("ensGene",genome,cachedir,sync=sync)
		en <- getUCSCTable("ensemblToGeneName",genome,cachedir,sync=sync)
		setnames(eg,"name","isoform.id")
		setnames(en,"name","isoform.id")
		setkey(eg,"isoform.id")
		setkey(en,"isoform.id")
		genes <- en[eg,list(chr=chrom,start=txStart+1,end=txEnd,strand=strand,name=value,gene.id=name2,isoform.id=isoform.id,cdsStart=cdsStart,cdsEnd=cdsEnd,exonCount=exonCount,exonStarts=exonStarts,exonEnds=exonEnds)]
		return(genes)
	}
}
# --------------------------------------------------------------------

# ====================================================================
# Internal Functions

# ====================================================================
