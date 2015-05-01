# Download, cache, and parse UCSC genome browser tables

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Load an annotation table from the UCSC Genome Browser as an R data.table
#'
#' If only \code{table} and \code{genome} are given, the function will load the data directly into the R workspace. If \code{cachedir} is a path to a directory, this directory will be used to maintain a cache of UCSC tables so they do not need to be re-downloaded on each call. If the data already exists and \code{sync=TRUE}, the function will only re-download and re-extract if the modified dates are different between the cachedir and remote copies. Note that start coordinates in these raw data tables are 0-based, whereas the Goldmine annotation functions will convert these to be 1-based.
#' @param table The UCSC string specific for the table to sync (e.g. "knownGene", "kgXref", etc)
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param cachedir A path to a directory where a cachedir cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly.
#' @param version If "latest" (default) then use the newest version of the table available. If set to a timestamp string of an archived table (format: YYYY-MM-DD-HH-MM-SS), then load this specific version. Obtain these strings by examining the file names under your cache directory. An archive file with a date stamp is saved automatically with each download of a new version. This feature only works if you have a cachedir cache that contains the desired versions.
#' @param sync If \code{TRUE}, then check if a newer version is available and download if it is. If \code{FALSE}, skip this check. Only has an effect if a cachedir cache directory (\code{cachedir}) is given.
#' @param url The root of the remote http URL to download UCSC data from (set by default to \code{http://hgdownload.cse.ucsc.edu/goldenPath/})
#; @param fread Use fread() from data.table and return a data.table object rather than the possibly much slower read.table() from base. Set to FALSE if you want a data.frame returned rather than a data.table. Default: TRUE
#' @return A data.frame or data.table of the desired UCSC table.
#' @export
getUCSCTable <- function(table, genome, cachedir=NULL, version="latest", sync=TRUE, url="http://hgdownload.cse.ucsc.edu/goldenPath/", fread=TRUE)
{
	# If we need to sync and a cachedir path has been given
	if((!is.null(cachedir))&(sync==TRUE))
	{
		# create dir (if we can) if it does not exist already
		if(!is.null(cachedir)){dir.create(cachedir, showWarnings = FALSE)}

		# check that our cachedir string has a trailing slash and add one if it does not
		cachedir <- normalizePath(cachedir)

		# check that the cachedir string points to a path that exists
		if(!file.exists(cachedir)){stop(paste("Error: Local path",cachedir,"does not exist. Please check the path and create the directory if necessary.",sep=" "))}
		syncUCSCTable(table, genome, url, cachedir)

		# Check if we are loading an archived version
		if(version!="latest")
		{
			cachedir.version=paste(".",version,sep="")
			print(paste("Using Archived Version: ",version,sep=" "))
		} else
		{
			cachedir.version=""
		}

		# Set paths
		cachedir.sql <- file.path(cachedir, genome, "database", paste(table, cachedir.version, ".sql", sep=""))
		cachedir.txt <- file.path(cachedir, genome, "database", paste(table, cachedir.version, ".txt", sep=""))
	} else if ((!is.null(cachedir))&(sync==FALSE))
	{
		# If we don't need to sync and a cachedir path has been given

		# Check if we are loading an archived version
		if(version!="latest")
		{
			cachedir.version=paste(".",version,sep="")
			print(paste("Using Archived Version: ",version,sep=" "))
		} else
		{
			cachedir.version=""
		}

		# Set paths
		cachedir.sql <- file.path(cachedir, genome, "database", paste(table, cachedir.version, ".sql", sep=""))
		cachedir.txt <- file.path(cachedir, genome, "database", paste(table, cachedir.version, ".txt", sep=""))

	} else if(is.null(cachedir))
	{
		# If cachedir path not given and we need to download to tempfile() on the fly

		url.dl.txt <- paste(url, genome , "/database/",table,".txt.gz",sep="")
		url.dl.sql <- paste(url, genome , "/database/",table,".sql",sep="")

		# Check that these URIs point to actual files (i.e. this table exists and we could download it if we want to)
		if((!RCurl::url.exists(url.dl.txt))|(!RCurl::url.exists(url.dl.sql))) {stop(paste("Error: Could not open table data URLs (",url.dl.txt," and ",url.dl.sql,"). Is the table name correct?",sep=""))}

		# Get tempfiles
		cachedir.txt.gz <- tempfile()
		cachedir.txt <- tempfile()
		cachedir.sql <- tempfile()

		# Download
		download.file(url.dl.txt, cachedir.txt.gz, quiet=FALSE)
		download.file(url.dl.sql, cachedir.sql, quiet=FALSE)

		# Gunzip the file
		R.utils::gunzip(cachedir.txt.gz, cachedir.txt, overwrite=TRUE, remove=TRUE)
	}

	# We can now open the data from the TXT
	if(fread==TRUE)
	{
		txt <- fread(cachedir.txt,sep="\t")
	} else
	{
		txt <- read.table(file=cachedir.txt, comment.char="", header=FALSE, stringsAsFactors=FALSE, sep="\t", quote="")
	}
	txt.nCols <- ncol(txt)

	# Parse SQL schema to get the row names
	mycols <- getTableHeaderFromSQL(cachedir.sql)

	# Set table with these row names
	if(fread==TRUE)
	{
		setnames(txt,mycols)
	} else
	{
		names(txt) <- mycols
	}

	# Return the dataframe
	txt
}
# --------------------------------------------------------------------

# ====================================================================

# ====================================================================
# Internal Functions

# --------------------------------------------------------------------
# internal function that does the syncing
syncUCSCTable <- function(table, genome, url, cachedir)
{
	# Generate URLs based on the given table name
	url.dl.txt <- paste(url, genome , "/database/",table,".txt.gz",sep="")
	url.dl.sql <- paste(url, genome , "/database/",table,".sql",sep="")

	# Check that these URIs point to actual files (i.e. this table exists and we could download it if we want to)
	if((!RCurl::url.exists(url.dl.txt))|(!RCurl::url.exists(url.dl.sql))) {stop(paste("Error: Could not open table data URLs (",url.dl.txt," and ",url.dl.sql,"). Is the table name correct?",sep=""))}

	# Create genome directory if it does not exist
	cachedir.dir.genome <- file.path(cachedir,genome)
	dir.create(cachedir.dir.genome, showWarnings = FALSE)

	# Create database directory if it does not exist
	cachedir.dir.database <- file.path(cachedir, genome, "database")
	dir.create(cachedir.dir.database, showWarnings = FALSE)

	# Generate full path of files to save cachedirly
	cachedir.file.txt <- file.path(cachedir, genome, "database", paste(table, ".txt.gz", sep=""))
	cachedir.file.sql <- file.path(cachedir, genome, "database", paste(table, ".sql", sep=""))
	cachedir.file.latest <- file.path(cachedir, genome, "database", paste(table, ".latest", sep=""))

	# ############################
	## Sync of table.txt.gz

	# Get version of last downloaded version
	if(file.exists(cachedir.file.latest))
	{
		file.latest <- readLines(cachedir.file.latest)
	} else
	{
		file.latest <- "file never downloaded"
	}

	# Get the remote file mod time from the HTTP header
	#mtime.remote <- httr::HEAD(url.dl.txt)$headers$`last-modified`
	#mtime.remote <- strptime(mtime.remote, format="%a, %d %b %Y %T", tz="GMT")
	# Noticed on 2014-12-23 that httr is mangling the modified header, this is not an ideal fix but restores function for now until httr can fix the bug
	mtime.remote <- toupper(httr::HEAD(url.dl.txt)$headers[3])
	#mtime.remote <- toupper(paste0(str_replace(names(mtime.remote)[3],"last-modified: ",""),":",mtime.remote[3]))
	file.latest <- toupper(file.latest)

	# Re-download only if this mod time is different than our own
	#if(is.na(file.mtime1.txt)){file.mtime1.txt="file never downloaded"}
	if(mtime.remote!=file.latest)
	{
		print(paste("Server version (", mtime.remote, ") is different than our latest (", file.latest,"), downloading new table...",sep=""))

		# Download the file
		download.file(url.dl.txt, cachedir.file.txt, quiet=FALSE)

		# Gunzip the file
		R.utils::gunzip(cachedir.file.txt, overwrite=TRUE, remove=TRUE)

		# I can't get any R downloaders to save the file with the server's modified date, so I am creating a file called table.latest.txt that stores the last modified string rather than working off what the filesystem reports
		# This will also help with an issue where MacOSX was touching these times and messing up the sync function!
		# The other advantage - we can delete the gzipped versions

		# Save datestamp of latest to a text file we can compare with later
		fc <- file(cachedir.file.latest)
		writeLines(mtime.remote, fc)
		close(fc)

		# Download table.sql
		# The timestamp is only checked for the .txt table, we pull a new SQL every time it changes to make things simpler (this way when the user wants a version the timestamps are the same between the SQL and the TXT)
		download.file(url.dl.sql, cachedir.file.sql, quiet=FALSE)

		# Create archive versions with datestamp in the file name of both the TXT and SQL files
		timestamp <- as.character(strptime(mtime.remote, format="%a, %d %b %Y %T", tz="GMT"))
		timestamp <- str_replace_all(str_replace_all(timestamp," ","-"),":","-")

		archive.file.txt <- file.path(cachedir, genome, "database", paste(table,timestamp, "txt", sep="."))
		archive.file.sql <- file.path(cachedir, genome, "database", paste(table,timestamp, "sql", sep="."))
		cachedir.file.txt2 <- file.path(cachedir, genome, "database", paste(table, ".txt", sep=""))

		file.copy(cachedir.file.txt2, archive.file.txt)
		file.copy(cachedir.file.sql, archive.file.sql)
		print(paste("Archived new TXT file to ", archive.file.txt, sep=""))
		print(paste("Archived new SQL file to ", archive.file.sql, sep=""))
	}

	# ############################
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# internal function to parse SQL to get table headers
getTableHeaderFromSQL <- function(sql.file)
{
	cols <- c()

	con  <- file(sql.file, open = "r")

	extract <- FALSE
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
	{
		linevec <- unlist(str_split(oneLine," "))
		#print(linevec[3])
		#if(is.na(linevec[3])){linevec[3]<-"KEY"}
		#print(linevec[3])

		if((extract==TRUE)&((linevec[3]=="KEY")|(linevec[3]=="PRIMARY")|(linevec[3]=="UNIQUE")|(linevec[3]=="DEFAULT")))
		{
			extract <- FALSE
		}
		if(extract==TRUE)
		{
			#print(linevec[3])
			colname <- linevec[3]
			colname <- str_replace_all(linevec[3],"`","")
			cols <- c(cols, colname)
		}
		if(linevec[1]=="CREATE")
		{
			extract <- TRUE
		}
	} 

	close(con)

	cols
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Table reading function
# Input: path to flatfile table
# Output: dataframe of table
read.table.ucsc <- function(path)
{
	#one function to read in all tables exported from UCSC the same way
	read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Table reading function for big tables - autodetects types from first 5 rows
# Input: path to flatfile table
# Output: dataframe of table
read.table.ucsc.big <- function(path)
{
	#detects colClasses on first 5 rows to speed up a big read in
	tab5rows <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", nrows = 5)
	classes <- sapply(tab5rows, class)
	tabAll <- read.table(file=path, comment.char="", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="", colClasses = classes)
	tabAll
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Read UCSC table files from disk and join all related tables
# Input: genome id, path with trailing "/"
# Output: joined annotation table dataframe (one row per gene isoform)
readUCSCAnnotation <- function(genome="hg19",path="")
{
	#load downloaded tables
	knownGene.path <- paste(path,genome,".knownGene.txt",sep="")
	#knownCanonical.path <- paste(path,genome,".knownCanonical.txt",sep="")
	kgXref.path <- paste(path,genome,".kgXref.txt",sep="")
	#cytoBand.path <- paste(path,genome,".cytoBand.txt",sep="")

	knownGene <- read.table.ucsc(knownGene.path)
	knownCanonical <- read.table.ucsc(knownCanonical.path)
	kgXref <- read.table.ucsc(kgXref.path)
	cytoBand <- read.table.ucsc.big(cytoBand.path)

	#join in gene symbols
	kg.sub <- knownGene
	names(kg.sub)[1] <- "name" 
	kgXref.sub <- kgXref
	names(kgXref.sub)[1] <- "name"

	kg.ann <- join(kg.sub,kgXref.sub,by="name",type="left")

	#join in canonical annotation (largest coding seq from nearest cluster)
	#kgCanon.sub <- data.frame(name=knownCanonical$transcript,canonical="1",stringsAsFactors=FALSE)
	
	#kg.ann.2 <- join(kg.ann.1,kgCanon.sub,by="name",type="left")

	#kg.ann <- kg.ann.2

	#join in cytological bands based on gene position
	#kg.ranges <- GRanges(seqnames=kg.ann$chrom,ranges=IRanges(start=kg.ann$txStart,end=kg.ann$txEnd),strand=kg.ann$strand,id=kg.ann$name)
	#cb.ranges <- with(cytoBand, GRanges(seqnames=X.chrom,ranges=IRanges(start=chromStart,end=chromEnd),band=name,gstain=gieStain))
	
	#band genes start in
	#kg.ranges$CytoBandStartIndex <- findOverlaps(kg.ranges,cb.ranges,select="first")
	#band genes end in (in case they span more than one)
	#kg.ranges$CytoBandEndIndex <- findOverlaps(kg.ranges,cb.ranges,select="last")
	#total bands spanned
	#kg.ranges$CytoBandsSpan <- countOverlaps(kg.ranges,cb.ranges)
	#create DF from GR metadata
	#cyto.map <- data.frame(name=kg.ranges$id,count=kg.ranges$CytoBandsSpan,start=kg.ranges$CytoBandStartIndex,end=kg.ranges$CytoBandEndIndex,stringsAsFactors=FALSE)

	#cyto.map$startname <- "NA"
	#cyto.map[is.na(cyto.map$start)==FALSE,]$startname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$start]$band
	#cyto.map$endname <- "NA"
	#cyto.map[is.na(cyto.map$start)==FALSE,]$endname <- cb.ranges[cyto.map[is.na(cyto.map$start)==FALSE,]$end]$band
	#cyto.map <- data.frame(name=cyto.map$name,cytoBandSpan=cyto.map$count,cytoBandStart=cyto.map$startname,cytoBandEnd=cyto.map$endname,stringsAsFactors=FALSE)

	#join cyto map data back into main annotation by id
	#kg.ann <- join(kg.ann,cyto.map,by="name",type="left")

	#add 1-based start coordinate column to remind viewer about this aspect of UCSC data
	#all ucsc start coords are 0 based and this code maintains that
	kg.ann$txStart.1based <- kg.ann$txStart + 1

	kg.ann$cdsStart.1based <- kg.ann$cdsStart + 1

	kg.ann
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Parse exon lists
# Input: data.frame from annotation reader
# Output: data.frame with each row representing an exon range
parseExons <- function(ann)
{
	#split by chrs to make it go faster
	chrs <- unique(ann$chrom)

	#my.ann <- ann[ann$chrom=="chr1",]	

	# parse out exon lists to give regions list where each individual exon is a range

	#ex <- foreach(j=1:length(chrs))
	ex <- foreach(j=1:length(chrs),.verbose=TRUE,.combine="rbind") %dopar%
	{
		my.ann <- ann[ann$chrom==chrs[j],]
		foreach(i=1:nrow(my.ann),.verbose=FALSE,.combine="rbind") %do%
		{
			chr <- my.ann[i,]$chrom
			starts <- as.numeric(unlist(strsplit(my.ann[i,]$exonStarts,",")))
			# correct for UCSC's 0-based system
			starts <- starts + 1
			ends <- as.numeric(unlist(strsplit(my.ann[i,]$exonEnds,",")))

			data.frame(chr=chr,start=starts,end=ends)
		}
	}


	ex
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
readRepeatMasker <- function(genome,path)
{
	rmsk.path <- paste(path,genome,".rmsk.txt",sep="")
	rmsk <- read.table.ucsc.big(rmsk.path)
	rmsk
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Counts of overlapping gene isoforms
# Input: genomic ranges object
# Output: vector of counts of how many gene isoforms overlap with region
getGenicOverlap <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# List of gene names overlapped by each range
# Input:
# Output: one column with a list of overlapping gene symbols for each row in the query regions
getGenicOverlapGenes <- function(regions.ranges, ann)
{
	ann.ranges <- with(ann, GRanges(seqnames=chrom,ranges=IRanges(start=txStart.1based,end=txEnd)))
	overlap <- findOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0

	overlap <- as.data.frame(overlap)
	overlap$name <- ann[overlap$subjectHits,]$geneSymbol

	out <- foreach(i=1:length(regions.ranges),.verbose=FALSE,.combine="c") %do%
	{
		hits <- overlap[overlap$queryHits==i,]
		genes <- unique(hits$name)
		paste(genes,collapse=", ")
	}

	out
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Find overlaps with upstream regions of genes
# Input: before = number bps upstream of TSS, after = number bps downstream of TSS
# Output:
getUpstreamOverlap <- function(regions.ranges, ann, before=1000, after=500)
{
	# add offsets, accounting for strandedness
	ups <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	ups$start.us <- NA
	ups$end.us <- NA

	ups[ups$strand=="+",]$start.us <- ups[ups$strand=="+",]$start - before
	ups[ups$strand=="+",]$end.us <- ups[ups$strand=="+",]$start + after


	ups[ups$strand=="-",]$start.us <- ups[ups$strand=="-",]$end - after
	ups[ups$strand=="-",]$end.us <- ups[ups$strand=="-",]$end + before

	ann.ranges <- with(ups, GRanges(seqnames=chr,ranges=IRanges(start=start.us,end=end.us)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Find overlaps with downstream regions of genes
# Input: before = number bps upstream of txEnd, after = number bps downstream of txEnd
# Output:
getDownstreamOverlap <- function(regions.ranges, ann, before=500, after=1000)
{
	# add offsets, accounting for strandedness
	downs <- with(ann,data.frame(chr=chrom,start=txStart.1based,end=txEnd,strand=strand))
	downs$start.ds <- NA
	downs$end.ds <- NA

	downs[downs$strand=="+",]$start.ds <- downs[downs$strand=="+",]$end - before
	downs[downs$strand=="+",]$end.ds <- downs[downs$strand=="+",]$end + after


	downs[downs$strand=="-",]$start.ds <- downs[downs$strand=="-",]$start - after
	downs[downs$strand=="-",]$end.ds <- downs[downs$strand=="-",]$start + before

	ann.ranges <- with(downs, GRanges(seqnames=chr,ranges=IRanges(start=start.ds,end=end.ds)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# 3' UTR Overlaps (gaps between cdsEnd and txEnd)
# Input:
# Output:
get3primeUTROverlap <- function(regions.ranges, ann)
{
	# filter for genes with a 3' UTR, accounting for strandedness
	utr <- with(ann,data.frame(chr=chrom, txStart=txStart.1based, txEnd=txEnd, strand=strand, cdsStart=cdsStart.1based, cdsEnd=cdsEnd))

	# filter non-coding transcripts which UCSC codes as cdsStart==cdsEnd
	utr <- utr[utr$cdsStart!=(utr$cdsEnd+1),]

	# filter out if cdsEnd == txEnd for (+) strand
	utr.p <- utr[(utr$strand=="+")&(utr$cdsEnd!=utr$txEnd),]

	# filter out if cdsStart == txStart for (-) strand because these are really the ends
	utr.m <- utr[(utr$strand=="-")&(utr$cdsStart!=utr$txStart),]

	# build list of utr regions
	utr.p$start.utr <- utr.p$cdsEnd
	utr.p$end.utr <- utr.p$txEnd

	utr.m$start.utr <- utr.m$txStart
	utr.m$end.utr <- utr.m$cdsStart

	reg <- rbind(utr.p,utr.m)

	ann.ranges <- with(reg, GRanges(seqnames=chr,ranges=IRanges(start=start.utr,end=end.utr)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Overlaps with exons
# Input: exon table from parseExons
# Output:
getExonOverlap <- function(regions.ranges, ann.ex)
{
	ann.ranges <- with(ann.ex, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
	overlap <- countOverlaps(regions.ranges,ann.ranges)
	#overlap[!is.na(overlap)] <- 1
	#overlap[is.na(overlap)] <- 0
	overlap
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getRepeatPercent <- function(regions.ranges,rmsk)
{
	#intersect with rmsk table, then collapse into nonoverlapping regions covering all repeats - want the % coverage of these regions versus the entire sequence length

	rmsk.ranges <- with(rmsk,GRanges(seqnames=genoName,ranges=IRanges(start=genoStart+1,end=genoEnd)))
	rmsk.ranges <- reduce(rmsk.ranges)

	#per <- foreach(i=1:length(regions.ranges),.combine="c") %do%
	per <- foreach(i=1:length(regions.ranges),.combine="c",.verbose=TRUE) %dopar%
	{
		print(i)
		cov <- intersect(regions.ranges[i],rmsk.ranges)
		#cov <- reduce(cov)
		#after reduction, the sum of the widths is the total coverage and we just need to divide this by the original width of the whole sequence
		rpt <- sum(width(cov)) / width(regions.ranges[i])
		rpt
	}

	per
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getRepeatPercentFast <- function(regions.ranges,genome,cachedir)
{
	#intersect with rmsk table, then collapse into nonoverlapping regions covering all repeats - want the % coverage of these regions versus the entire sequence length
	
	# Get from UCSC table
	rmsk <- getUCSCTable("rmsk",genome,cachedir)

	# Make into GRanges
	rmsk.ranges <- with(rmsk,GRanges(seqnames=genoName,ranges=IRanges(start=genoStart+1,end=genoEnd)))
	rmsk.ranges <- reduce(rmsk.ranges)

	per <- calcPercentOverlap(regions.ranges,rmsk.ranges)

	per
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getDistTSS <- function(regions.ranges,ann)
{
	# create ranges object with just the TSS
	# need to use txEnd for genes on the "-" strand
	ann.starts <- with(ann,data.frame(chrom,txStart.1based,txEnd,strand))

	ann.starts$tss <- NA
	ann.starts[ann.starts$strand=="+",]$tss <- ann.starts[ann.starts$strand=="+",]$txStart.1based
	ann.starts[ann.starts$strand=="-",]$tss <- ann.starts[ann.starts$strand=="-",]$txEnd

	ann.ranges <- with(ann.starts, GRanges(seqnames=chrom,ranges=IRanges(start=tss,end=tss)))

	#overlap <- countOverlaps(regions.ranges,ann.ranges)

	# if the overlaps number is zero, we need to find the nearest region

	#nearest(regions.ranges,ann.ranges)

	dtss <- as.data.frame(distanceToNearest(regions.ranges,ann.ranges))
	dtss[,3]
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#calculates from the center of the DMS rather than the nearest outer bound
getDistTSSCenter <- function(regions.ranges, genome, cachedir)
{
	# create ranges object with just the TSS
	# need to use txEnd for genes on the "-" strand
	#ann.starts <- with(ann,data.frame(chrom,txStart.1based,txEnd,strand))
	ann.starts <- getUCSCTable("knownGene",genome,cachedir)
	ann.starts[strand=="+",tss:=txStart+1]
	ann.starts[strand=="-",tss:=txEnd]

	ann.ranges <- with(ann.starts, GRanges(seqnames=chrom,ranges=IRanges(start=tss,end=tss)))

	#overlap <- countOverlaps(regions.ranges,ann.ranges)

	# if the overlaps number is zero, we need to find the nearest region

	#nearest(regions.ranges,ann.ranges)

	centers <- round(start(regions.ranges)+((end(regions.ranges) - start(regions.ranges))/2),digits=0)

	centers.ranges <- GRanges(seqnames=seqnames(regions.ranges),ranges=IRanges(start=centers,end=centers))

	dtss <- as.data.frame(distanceToNearest(centers.ranges,ann.ranges))
	dtss[,3]
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getDistTSE <- function(regions.ranges,genome,cachedir)
{
	# create ranges object with just the TSS
	# need to use txEnd for genes on the "-" strand
	ann.starts <- with(ann,data.frame(chrom,txStart.1based,txEnd,strand))

	ann.starts$tse <- NA
	ann.starts[ann.starts$strand=="+",]$tse <- ann.starts[ann.starts$strand=="+",]$txEnd
	ann.starts[ann.starts$strand=="-",]$tse <- ann.starts[ann.starts$strand=="-",]$txStart.1based

	ann.ranges <- with(ann.starts, GRanges(seqnames=chrom,ranges=IRanges(start=tse,end=tse)))

	#overlap <- countOverlaps(regions.ranges,ann.ranges)

	# if the overlaps number is zero, we need to find the nearest region

	#nearest(regions.ranges,ann.ranges)

	dtse <- as.data.frame(distanceToNearest(regions.ranges,ann.ranges))
	dtse[,3]
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getDistTSECenter <- function(regions.ranges,genome,cachedir)
{
	# create ranges object with just the TSS
	# need to use txEnd for genes on the "-" strand
	ann.starts <- getUCSCTable("knownGene",genome,cachedir)
	ann.starts[strand=="-",tse:=txStart+1]
	ann.starts[strand=="+",tse:=txEnd]
	ann.ranges <- with(ann.starts, GRanges(seqnames=chrom,ranges=IRanges(start=tse,end=tse)))

	#overlap <- countOverlaps(regions.ranges,ann.ranges)

	# if the overlaps number is zero, we need to find the nearest region

	#nearest(regions.ranges,ann.ranges)

	centers <- round(start(regions.ranges)+((end(regions.ranges) - start(regions.ranges))/2),digits=0)

	centers.ranges <- GRanges(seqnames=seqnames(regions.ranges),ranges=IRanges(start=centers,end=centers))


	dtse <- as.data.frame(distanceToNearest(centers.ranges,ann.ranges))
	dtse[,3]
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
getFreqCpG <- function(seq)
{
	dinucleotideFrequency(seq,as.prob=TRUE)[,7]
}
# --------------------------------------------------------------------

# ====================================================================



