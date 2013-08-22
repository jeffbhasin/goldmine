\name{annotateFeatures}
\alias{annotateFeatures}
\title{Annotate percent of query that overlaps with any regions in features, or report unique IDs from features for a specified column that overlap with query}
\usage{
  annotateFeatures(query, features, getid = NULL,
    ucsc = FALSE)
}
\arguments{
  \item{query}{A data.frame of regions to annotate. Must
  contain the columns "chr", "start", "end", and the
  "start" coordinates must be 1-based. Columns will not be
  retained in output, but the output vector will match
  query row for row.}

  \item{features}{A data.frame of regions to overlap query
  against. Must contain the columns "chr", "start", "end",
  and the "start" coordinates must be 1-based. Or, must be
  a raw UCSC table from getUCSCTable() if ucsc=TRUE is also
  given.}

  \item{getid}{If NULL, return percentage overlaps. If a
  string, the name of a column in features to pull
  overlapping IDs or other data fields from.}

  \item{ucsc}{Set to TRUE if the table is a raw UCSC table
  fetched with getUCSCTable()}
}
\value{
  A vector of the requested overlap data, matching row for
  row with query.
}
\description{
  Given a set of query regions in a data.frame containing
  the columns "chr", "start", and "end", returns a vector
  in the same length and order as the rows in query with
  either the percent of each query that overlap with any
  regions in features (default), or, if getid is given as a
  string, the function will report a comma separated string
  of IDs that are involved. Raw UCSC tables can be passed
  directly if ucsc=TRUE is given, and this option will
  adjust for 0-based start coordinates and UCSC's column
  names for the coordinates.
}

