# goldmine: Genome context annotation in R or from the command line

## About
* Obtains data by direct downloading and updating of a local mirror of select UCSC Genome Browser annotation tables
* R package contains functions to assess genomic context of any given set of genomic ranges by performing overlaps with regions of annotated genomic features and produce long, short, and plot outputs
* Also contains a script that allows annotation of a flat file from the command line, without needing the user to code in R

## Development Progress
* getUCSCTable has been implemented which supports downloading and loading of any UCSC table from any UCSC genome into R, and can optionally sync a local cache of these tables
* goldmine() has been implemented which creates a list of tables describing genomic context with respect to genes and features
* Functions for differential testing are being updated to integrate with the goldmine() output list

## Installation
From R:
```r
library(devtools)
install_github("gnomegenome/goldmine")
```
## Usage
Example: Loading knownGene on the fly to a data.frame
```r
library(goldmine)
knownGene <- getUCSCTable("knownGene","hg19")
```
