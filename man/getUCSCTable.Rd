\name{getUCSCTable}
\alias{getUCSCTable}
\title{Load an annotation table from the UCSC Genome Browser as an R data.frame}
\usage{
  getUCSCTable(table, genome, local = NULL,
    version = "latest", sync = TRUE,
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/")
}
\arguments{
  \item{table}{The UCSC string specific for the table to
  sync (e.g. "knownGene", "kgXref", etc)}

  \item{genome}{The UCSC string specific to the genome to
  be downloaded (e.g. "hg19", "hg19", "mm10", etc)}

  \item{local}{A path to a directory where a local cache of
  UCSC tables are stored. If equal to \code{NULL}
  (default), the data will be downloaded to temporary files
  and loaded on the fly.}

  \item{version}{If "latest" (default) then use the newest
  version of the table available. If set to a timestamp
  string of an archived table (format:
  YYYY-MM-DD-HH-MM-SS), then load this specific version.
  Obtain these strings by examining the file names under
  your cache directory. An archive file with a date stamp
  is saved automatically with each download of a new
  version. This feature only works if you have a local
  cache that contains the desired versions.}

  \item{sync}{If \code{TRUE}, then check if a newer version
  is available and download if it is. If \code{FALSE}, skip
  this check. Only has an effect if a local cache directory
  (\code{local}) is given.}

  \item{url}{The root of the remote http URL to download
  UCSC data from (set by default to
  \code{http://hgdownload.cse.ucsc.edu/goldenPath/})}
}
\value{
  A data.frame of the desired UCSC table.
}
\description{
  If only \code{table} and \code{genome} are given, the
  function will load the data directly into the R
  workspace. If \code{local} is a path to a directory, this
  directory will be used to maintain a local cache of UCSC
  tables so they do not need to be re-downloaded on each
  call. If the data already exists and \code{sync=TRUE},
  the function will only re-download and re-extract if the
  modified dates are different between the local and remote
  copies.
}

