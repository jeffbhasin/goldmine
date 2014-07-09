# #############################################################################
# Download, update, and parse UCSC genome browser tables
# Author: Jeffrey Bhasin <jeffb@case.edu>
# Created: 2013-07-31
# #############################################################################

# =============================================================================
# User-Facing Functions

# -----------------------------------------------------------------------------
#' Load an annotation table from the UCSC Genome Browser as an R data.frame
#'
#' If only \code{table} and \code{genome} are given, the function will load the data directly into the R workspace. If \code{cachedir} is a path to a directory, this directory will be used to maintain a cachedir cache of UCSC tables so they do not need to be re-downloaded on each call. If the data already exists and \code{sync=TRUE}, the function will only re-download and re-extract if the modified dates are different between the cachedir and remote copies.
#' @param table The UCSC string specific for the table to sync (e.g. "knownGene", "kgXref", etc)
#' @param genome The UCSC string specific to the genome to be downloaded (e.g. "hg19", "hg19", "mm10", etc)
#' @param cachedir A path to a directory where a cachedir cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly.
#' @param version If "latest" (default) then use the newest version of the table available. If set to a timestamp string of an archived table (format: YYYY-MM-DD-HH-MM-SS), then load this specific version. Obtain these strings by examining the file names under your cache directory. An archive file with a date stamp is saved automatically with each download of a new version. This feature only works if you have a cachedir cache that contains the desired versions.
#' @param sync If \code{TRUE}, then check if a newer version is available and download if it is. If \code{FALSE}, skip this check. Only has an effect if a cachedir cache directory (\code{cachedir}) is given.
#' @param url The root of the remote http URL to download UCSC data from (set by default to \code{http://hgdownload.cse.ucsc.edu/goldenPath/})
#; @param fread Use fread() from data.table and return a data.table object rather than the possibly much slower read.table() from base. Set to FALSE if you want a data.frame returned rather than a data.table. Default: TRUE
#' @return A data.frame of the desired UCSC table.
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
		if((!url.exists(url.dl.txt))|(!url.exists(url.dl.sql))) {stop(paste("Error: Could not open table data URLs (",url.dl.txt," and ",url.dl.sql,"). Is the table name correct?",sep=""))}

		# Get tempfiles
		cachedir.txt.gz <- tempfile()
		cachedir.txt <- tempfile()
		cachedir.sql <- tempfile()

		# Download
		download.file(url.dl.txt, cachedir.txt.gz, quiet=FALSE)
		download.file(url.dl.sql, cachedir.sql, quiet=FALSE)

		# Gunzip the file
		gunzip(cachedir.txt.gz, cachedir.txt, overwrite=TRUE, remove=TRUE)
	}

	# We can now open the data from the TXT
	if(fread==TRUE)
	{
		txt <- fread(cachedir.txt)
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
# -----------------------------------------------------------------------------
# =============================================================================

# =============================================================================
# Internal Functions


# =============================================================================

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
	mtime.remote <- httr::HEAD(url.dl.txt)$headers$`last-modified`
	#mtime.remote <- strptime(mtime.remote, format="%a, %d %b %Y %T", tz="GMT")

	# Re-download only if this mod time is different than our own
	#if(is.na(file.mtime1.txt)){file.mtime1.txt="file never downloaded"}
	if(mtime.remote!=file.latest)
	{
		print(paste("Server version (", mtime.remote, ") is different than our latest (", file.latest,"), downloading new table...",sep=""))

		# Download the file
		download.file(url.dl.txt, cachedir.file.txt, quiet=FALSE)

		# Gunzip the file
		gunzip(cachedir.file.txt, overwrite=TRUE, remove=TRUE)

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