# Annotate genomic ranges using reference data about genes and features

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Explore relationships between a set of genomic ranges and known genes/features
#'
#' Computes the overlap between a query set of genomic ranges given as a GenomicRanges, data.frame, or data.table with gene and feature sets of interest. Reports both summarized overlaps (same number of rows as the query - a "wide format") and in separate tables, individual overlap events (one row for each pair of overlapping query and gene/feature item - a "long format" similar to an inner join).
#' @param query A GenomicRanges, data.frame, or data.table of regions to annotate. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based. All additional columns will be retained in the output object.
#' @param genes Genes of interest from the output table of getGenes(). If not given, will default to the UCSC knownGene table.
#' @param features A list() of GenomicRanges, data.table, or data.frame objects giving feature sets of interest. If a data.frame or data.table, must contain the columns "chr", "start", "end", where the "start" coordinates are 1-based. All additional columns will be retained in the output object. See also the getFeatures() function.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param promoter A numeric vector of length 2 specifying the number of bp upstream and downstream of transcription start sites for which to create promoter ranges. Given as c(upstream,downstream). Note that "upstream" in the context of the 5' end of the gene means out from the gene body.
#' @param end3 A numeric vector of length 2 specifying the number of bp upstream and downstream of transcription end sites for which to create gene 3' end ranges. Given as c(upstream,downstream). Note that "upstream" in the context of the 3' end of the gene means into the gene body.
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @param sync If TRUE, then check if newer versions of UCSC tables are available and download them if so. If FALSE, skip this check. Can be used to freeze data versions in an analysis-specific cachedir for reproducibility.
#' @return A list: "context" shows a percent overlap for each range in the query set with gene model regions and each feature set ("wide" format - same number of rows as the query and in the same order), "genes" contains a detailed view of each query region overlap with individual gene isoforms ("long" format - a row for each pair of query and isoform overlaps), "features" is a list of tables which for each table given in the "features" argument which contain a row for each instance of a query region overlapping with a feature region (also "long" format).
#' @export
goldmine <- function(query, genes=getGenes(geneset="ucsc", genome=genome, cachedir=cachedir), features=list(), promoter=c(1000,500), end3=c(1000,1000), contextonly=FALSE, genome, cachedir, sync=TRUE)
{
	# Validate and convert query input
	query.gr <- makeGRanges(query)
	#query.dt <- makeDT(query)
	#genes.dt <- makeDT(genes)
	genes.gr <- makeGRanges(genes, strand=T)
	if(!is.list(features)){stop("features must be given as a list. set as list() to disable feature overlaps.")}
	features.gr <- lapply(features,makeGRanges)
	#features.dt <- lapply(features,makeDT)

	# Get chromosome lengths
	chromInfo <- getUCSCTable("chromInfo", genome, cachedir, sync=sync)
	chromInfo$chr <- chromInfo$chrom

	# Set key ID fields for easy joins
	#genes.dt[,srow:=1:nrow(genes.dt)]
	#query.dt[,qrow:=1:nrow(query.dt)]
	genes.gr$srow <- 1:length(genes.gr)
	query.gr$qrow <- 1:length(query.gr)

	# Extract gene models
	genemodels <- suppressWarnings(goldmine:::getGeneModels(genes=genes.gr, promoter=promoter, end3=end3, genome=genome, cachedir=cachedir, sync=sync))

	# Do the context annotation ("wide format" - returns same rows as original plus annotation columns)
	message("Generating context annotation - genes")
	genemodels.per <- suppressWarnings(lapply(genemodels,function(x) goldmine:::calcPercentOverlap(query.gr,reduce(x))))
	names(genemodels.per) <- paste0(names(genemodels.per),"_per")
	ann.gene <- do.call(cbind,genemodels.per)
	ann <- addNearest(query.gr,genes.gr,id="name",prefix="genes")
	ann <- cbind(ann, ann.gene)

	# Call categories (mutually exclusive)
	p <- ann$promoter_per > 0
	i <- ann$intron_per > 0
	e <- ann$exon_per > 0
	t <- ann$end3_per > 0
	calls <- rep("intergenic",nrow(ann))
	calls[i] <- "intron"
	calls[e] <- "exon"
	calls[t] <- "3' end"
	calls[p] <- "promoter"
	ann$call <- calls

	# Want two gene columns - one that tells the gene(s) associated with the context calls and one that gives all overlapping genes OR the nearest gene
	# If intergenic, just leave as nearest gene
	ann$call_genes <- "" 
	ann$overlapped_genes <- ""

	genemodels$end3$num <- 0
	genemodels$promoter$num <- 0
	genemodels$end3$con <- "3' end"
	genemodels$promoter$con <- "promoter"
	genemodels$intron$con <- "intron"
	genemodels$exon$con <- "exon"
	gparts <- suppressWarnings(c(genemodels$exon,genemodels$intron,genemodels$promoter,genemodels$end3))
	
	# call_genes are the genes behind the context call
	fo <- as.data.frame(findOverlaps(query.gr,gparts))
	fo$call <- ann[fo$queryHits,]$call
	fo$gpart <- gparts[fo$subjectHits]$con
	fo$gene <- genes.gr[gparts[fo$subjectHits]$srow]$name
	focall <- fo[fo$call==fo$gpart,]
	focall <- data.table(focall)
	focall <- focall[,list(call_genes=toString(unique(gene))),by="queryHits"]
	ann[focall$queryHits,]$call_genes <- focall$call_genes

	# associated_genes are nearest genes for intergenic, and all genes for which any part is overlapped for the others
	fo <- data.table(fo)
	fo <- fo[,list(associated_genes=toString(unique(gene))),by="queryHits"]

	ann[fo$queryHits,]$overlapped_genes <- fo$associated_genes

	# order all the gene columns together
	g <- ann$genes_nearest
	d <- ann$genes_dist
	ann[,genes_nearest:=NULL]
	ann[,genes_dist:=NULL]

	ann$nearest_genes <- g
	ann$distance_to_nearest_gene <- d

	if(length(features)>0)
	{
		message("Generating context annotation - features")
		features.per <- lapply(features.gr,function(x) goldmine:::calcPercentOverlap(query.gr,reduce(x)))
		names(features.per) <- paste0(names(features.per),"_per")
		ann <- cbind(ann, do.call(cbind,features.per))
	}

	# Add URL
	ann$url <- goldmine:::getBrowserURLs(query.gr,genome)

	if(contextonly==FALSE)
	{
			# Do the gene overlaps if asked for
			message("Generating genes report")
			reportGenes2 <- function(query.gr, genes.gr)
			{
				chrs <- as.character(seqnames(query.gr))
				starts <- start(query.gr)
				ends <- end(query.gr)

				# Initial Intersect
				fo <- findOverlaps(query.gr, genes.gr)
				foq <- queryHits(fo)
				fos <- subjectHits(fo)

				# Make initial joined columns
				qrow=query.gr[foq]$qrow
				srow=genes.gr[fos]$srow
				query.chr=chrs[foq]
				query.start=starts[foq]
				query.end=ends[foq]
				gene.symbol=genes.gr[fos]$name
				gene.id=genes.gr[fos]$gene.id
				isoform.id=genes.gr[fos]$isoform.id
				isoform.chr=as.character(seqnames(genes.gr[fos]))
				isoform.start=start(genes.gr[fos])
				isoform.end=end(genes.gr[fos])
				isoform.strand=as.character(strand(genes.gr[fos]))
				overlap.bp=suppressWarnings(goldmine:::calcPercentOverlap(query.gr, genes.gr, sum.all=FALSE, report.bp=TRUE))
				query.overlap.per=suppressWarnings(goldmine:::calcPercentOverlap(query.gr, genes.gr, sum.all=FALSE))
				isoform.overlap.per=suppressWarnings(goldmine:::calcPercentOverlap(genes.gr, query.gr, sum.all=FALSE))

				anng <- data.table(qrow=qrow, srow=srow, query.chr=query.chr, query.start=query.start, query.end=query.end, gene.symbol=gene.symbol, gene.id=gene.id, isoform.id=isoform.id, isoform.chr=isoform.chr, isoform.start=isoform.start, isoform.end=isoform.end, isoform.strand=isoform.strand, overlap.bp=overlap.bp, query.overlap.per=query.overlap.per, isoform.overlap.per=isoform.overlap.per)

				# add noncoding flag
				anng$noncoding <- genes.gr[fos]$cdsStart==genes.gr[fos]$cdsEnd

				# Make isoform-specific percent overlap columns
				# Promoter
				prom.o <- goldmine:::calcOverlapForReport(query.gr,genemodels$promoter, all=FALSE)
				colnames(prom.o) <- c("qrow", "srow", "Promoter")
				setkeyv(anng,c("qrow","srow"))
				prom.o <- data.table(prom.o)
				setkeyv(prom.o,c("qrow","srow"))
				anng <- merge(anng,prom.o,all.x=T,all.y=F)
				anng[is.na(anng$"Promoter"),Promoter:=0]

				# Exon and Intron Diagram
				exon.o <- goldmine:::calcOverlapForReportExons(query.gr, genemodels$exon)
				intron.o <- goldmine:::calcOverlapForReportExons(query.gr, genemodels$intron)
				if(nrow(exon.o)>0)
				{
					exon.o$type <- "E"
				}
				if(nrow(intron.o)>0)
				{
					intron.o$type <- "I"
				}
				all.o <- rbind(exon.o, intron.o)

				message("Generating exon/intron overlap diagrams")
				kgxs <- rbind(data.table(srow=genemodels$exon$srow, num=genemodels$exon$num, type="E"), data.table(srow=genemodels$intron$srow, num=genemodels$intron$num, type="I"))
				pairs.dt <- data.table(qrow=anng$qrow, srow=anng$srow)
				setkey(pairs.dt, qrow, srow)
				kgxs.dt <- data.table(kgxs)
				setkey(kgxs.dt, srow)
				kgx.dt <- merge(pairs.dt, kgxs.dt, by="srow", allow.cartesian=TRUE)

				# Join our percents in
				all.o.dt <- data.table(all.o)
				setkey(all.o.dt, qrow, srow, num, type)
				ts.dt <- merge(kgx.dt, all.o.dt, by=c("qrow","srow","num","type"), all.x=TRUE, all.y=FALSE)

				# Assign the rest as 0%
				if(sum(is.na(ts.dt$per))>0)
				{
					ts.dt[is.na(ts.dt$per),per:=0]
				}
				ts.dt[,strand:=as.character(strand(genes.gr))[srow]]

				# Reverse numbers if strand="-"
				ts1 <- ts.dt[strand=="-",list(num=rev(num),type=type,per=per,strand=strand),by=c("qrow","srow")]
				ts2 <- ts.dt[strand=="+",]
				ts.dt <- rbind(ts1,ts2)

				ts.dt <- ts.dt[,list(qrow=qrow, srow=srow, strand=strand, string=paste0(type,num," (",per,")"))]

				# Reverse the toString if the strand is "-"

				#ts.dt <- ts.dt[,list("ExonIntron"=toString(string)),by=c("qrow","srow")]
				st.dt1 <- ts.dt[strand=="+",list(ExonIntron=toString(string)),by=c("qrow","srow")]
				st.dt2 <- ts.dt[strand=="-",list(ExonIntron=toString(rev(string))),by=c("qrow","srow")]
				st.dt <- rbind(st.dt1,st.dt2)

				anng <- merge(anng,st.dt,by=c("qrow","srow"),all.x=T)

				# 3' Ends
				end.o <- goldmine:::calcOverlapForReport(query.gr,genemodels$end3, all=FALSE)
				colnames(end.o) <- c("qrow", "srow", "3' End")
				setkeyv(anng,c("qrow","srow"))
				end.o <- data.table(end.o)
				setkeyv(end.o,c("qrow","srow"))
				anng <- merge(anng,end.o,all.x=T,all.y=F)
				anng[is.na(anng$"3' End"),"3' End":=0]

				# add browser links if this was requested
				anng$url <- goldmine:::getBrowserURLs(GRanges(seqnames=anng$query.chr, ranges=IRanges(start=anng$query.start,end=anng$query.end)), genome)

				return(anng)

			}
			rg <- reportGenes2(query.gr, genes.gr)

			# Do the feature overlaps if asked for, separate for each feature on the list
			if(length(features)>0)
			{
			message("Generating features report")
			reportFeatures2 <- function(query.gr, x.gr)
			{
				features.chr <- as.character(seqnames(x.gr))
				features.start <- start(x.gr)
				features.end <- end(x.gr)
				chrs <- as.character(seqnames(query.gr))
				starts <- start(query.gr)
				ends <- end(query.gr)
				query.gr$qrow <- 1:length(query.gr)
				x.gr$srow <- 1:length(x.gr)
				overs <- goldmine:::calcOverlapForReport(query.gr, x.gr)

				annf <- data.table(query.chr=chrs[overs$qrow], query.start=starts[overs$qrow], query.end=ends[overs$qrow], feature.chr=features.chr[overs$srow], feature.start=features.start[overs$srow], feature.end=features.end[overs$srow], overlap.query.per=overs$per, overlap.feature.per=overs$sper, overlap.bp=overs$bp)
				annq <- as.data.frame(values(query.gr)[overs$qrow,])
				colnames(annq) <- paste("query",colnames(annq),sep="_")
				anns <- as.data.frame(values(x.gr)[overs$srow,])
				colnames(anns) <- paste("feature",colnames(anns),sep="_")

				annf2 <- cbind(annf, annq, anns)
				return(annf2)
			}
			rf <- lapply(features.gr,function(x) reportFeatures2(query.gr, x))
			} else
			{
				rf <- list()
			}
			return(list(context=ann, genes=rg, features=rf))
	} else
	{
		return(ann)
	}
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Write individual CSV files to disk from the output of goldmine()
#'
#' Write a CSV file for each output table in a goldmine() output list object. Useful for importing to spreadsheet applications. Will save the "context", "genes", and any tables in "features" to individual CSV files.
#' @param gm The output list object from goldmine().
#' @param path The directory to write the files into (default: current working directory).
#' @export
gmWrite <- function(gm,path=".")
{
	if(!file.exists(path)){dir.create(path)}

	f <- paste0(path,"/context.csv")
	message("Writing CSV File: ",f)
	write.csv(gm$context,file=f,row.names=F)
	f <- paste0(path,"/genes.csv")
	message("Writing CSV File: ",f)
	write.csv(gm$genes,file=f,row.names=F)
	if(length(gm$features)>0)
	{
		for(i in 1:length(gm$features))
		{
			f <- paste0(path,"/features_",names(gm$features)[i],".csv")
			message("Writing CSV File: ",f)
			write.csv(gm$features[[i]],file=f,row.names=F)
		}
	}
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Return sets of ranges for individual gene model components based on the contents of the output from getGenes()
#'
#' For custom analysis that requires the genomic ranges of gene model components.
#' @param genes Genes of interest from the output table of getGenes(). If not given, will default to the UCSC knownGene table.
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param promoter A numeric vector of length 2 specifying the number of bp upstream and downstream of transcription start sites for which to create promoter ranges. Given as c(upstream,downstream). Note that "upstream" in the context of the 5' end of the gene means out from the gene body.
#' @param end3 A numeric vector of length 2 specifying the number of bp upstream and downstream of transcription end sites for which to create gene 3' end ranges. Given as c(upstream,downstream). Note that "upstream" in the context of the 3' end of the gene means into the gene body.
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @param sync If TRUE, then check if newer versions of UCSC tables are available and download them if so. If FALSE, skip this check. Can be used to freeze data versions in an analysis-specific cachedir for reproducibility.
#' @return A list containing one GenomicRanges object for each of the gene model portions: Promoters, 3' Ends, Exons, Introns, Intergenic, 5' UTRs, 3' UTRs. The "srow" column can be used to match individual ranges to individual genes in the table given to the "genes" argument by row number.
#' @export
getGeneModels <- function(genes=getGenes(geneset="ucsc", genome=genome, cachedir=cachedir), promoter=c(1000,500), end3=c(1000,1000), genome, cachedir, sync=TRUE)
{
	message("Computing gene models")

	genes.gr <- makeGRanges(genes,strand=T)
	genes.gr$srow <- 1:length(genes.gr)
	chromInfo <- getUCSCTable("chromInfo", genome, cachedir, sync=sync)
	chromInfo$chr <- chromInfo$chrom
	seqlengths(genes.gr) <- chromInfo[match(seqlevels(genes.gr), chromInfo$chrom),]$size

	genemodels <- list()

	# Promoters
	prom.gr <- suppressWarnings(promoters(genes.gr,upstream=promoter[1],downstream=promoter[2]))
	dat <- values(prom.gr)$srow
	values(prom.gr) <- NULL
	prom.gr$srow <- dat
	genemodels$promoter <- prom.gr

	# 3' Ends
	# Need to flip strand so now "promoters" works from the end of the gene
	tmp.gr <- genes.gr
	strand(tmp.gr)[strand(genes.gr)=="+"] <- "-"
	strand(tmp.gr)[strand(genes.gr)=="-"] <- "+"
	ends.gr <- suppressWarnings(promoters(tmp.gr,upstream=end3[2],downstream=end3[1]))
	dat <- values(ends.gr)$srow
	values(ends.gr) <- NULL
	ends.gr$srow <- dat
	strand(ends.gr) <- strand(genes.gr)
	genemodels$end3 <- ends.gr

	# Exons
	exonstarts <- str_replace(genes.gr$exonStarts,",$","")
	exonstartvec <- unlist(str_split(exonstarts,","))
	exonstartvec <- as.numeric(exonstartvec)+1
	exonends <- str_replace(genes.gr$exonEnds,",$","")
	exonendvec <- as.numeric(unlist(str_split(exonends,",")))
	exonchrs <- rep(as.character(seqnames(genes.gr)), genes.gr$exonCount)
	exonsrows <- rep(genes.gr$srow, genes.gr$exonCount)
	exonnums <- unlist(sapply(genes.gr$exonCount, FUN=function(x) seq(1, x)))

	exon.gr <- GRanges(seqnames=exonchrs, ranges=IRanges(start=exonstartvec,end=exonendvec), srow=exonsrows, num=exonnums)
	genemodels$exon <- exon.gr

	# Introns
	exontrans <- rep(genes.gr$srow, genes.gr$exonCount)
	int.dt <- data.table(chr=exonchrs[-1],startiso=exontrans[-1],exonstart=exonstartvec[-1],endiso=exontrans[1:(length(exontrans)-1)],exonend=exonendvec[1:(length(exonendvec)-1)],exonnums=exonnums[-1]-1)
	tron2 <- int.dt[startiso==endiso,]
	intron.gr <- GRanges(tron2$chr,IRanges(tron2$exonend,tron2$exonstart-1), srow=tron2$startiso, num=tron2$exonnums)
	genemodels$intron <- intron.gr

	# Intergenic
	intergenic.gr <- GRanges(seqnames=chromInfo$chrom, ranges=IRanges(start=1, end=chromInfo$size))
	seqlengths(intergenic.gr) <- chromInfo[match(names(seqlengths(intergenic.gr)), chromInfo$chr),]$size
	intergenic.gr <- setdiff(intergenic.gr, prom.gr)
	intergenic.gr <- setdiff(intergenic.gr, ends.gr)
	intergenic.gr <- setdiff(intergenic.gr, exon.gr)
	intergenic.gr <- setdiff(intergenic.gr, intron.gr)
	genemodels$intergenic <- intergenic.gr

	# Also report UTRs if asked for?
	# Make 5' UTR
	kg <- makeDT(genes)
	# distance between txStart and cdsStart for + genes
	# distance between txEnd and cdsEnd for - genes
	kg.p <- kg[(kg$strand=="+")&(kg$start!=kg$cdsStart)&(kg$cdsEnd!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$end!=kg$cdsEnd)&(kg$cdsEnd!=kg$cdsStart),]

	utr5.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chr, ranges=IRanges(start=start, end=cdsStart))),with(kg.m, GRanges(seqnames=kg.m$chr, ranges=IRanges(start=cdsEnd, end=end)))))

	# Make 3' UTR
	kg.p <- kg[(kg$strand=="+")&(kg$end!=kg$cdsEnd)&(kg$cdsEnd!=kg$cdsStart),]
	kg.m <- kg[(kg$strand=="-")&(kg$start!=kg$cdsStart)&(kg$cdsEnd!=kg$cdsStart),]

	utr3.gr <- suppressWarnings(c(with(kg.p, GRanges(seqnames=kg.p$chr, ranges=IRanges(start=cdsEnd, end=end))),with(kg.m, GRanges(seqnames=kg.m$chr, ranges=IRanges(start=start, end=cdsStart)))))

	genemodels$utr5 <- utr5.gr
	genemodels$utr3 <- utr3.gr

	return(genemodels)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Obtain feature sets from UCSC genome browser tables
#'
#' Given a vector of table names from the UCSC genome browser that all contain "chrom", "chromStart", and "chromEnd" fields, converts them to input suitable for the goldmine() "features" argument (changes column names and adjusts 0-based start coordinates to 1-based).
#' @param tables A vector of table names from UCSC (default: DnaseI and TFBS from encode for hg19).
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @param sync If TRUE, then check if newer versions of UCSC tables are available and download them if so. If FALSE, skip this check. Can be used to freeze data versions in an analysis-specific cachedir for reproducibility.
#' @export
getFeatures <- function(tables=c("wgEncodeRegDnaseClusteredV3","wgEncodeRegTfbsClusteredV3"), genome, cachedir, sync=TRUE)
{
	#if(genome!="hg19"){stop("This shortcut function is designed for genome hg19 only. Pleasure use getUCSCTable() to build feature sets for any other UCSC genome as desired.")}
	#tables <- c("wgEncodeRegDnaseClusteredV2","wgEncodeRegTfbsClusteredV3","tfbsConsSites","gwasCatalog", "pubsBlat", "cosmic", "oreganno", "vistaEnhancers", "phastConsElements100way")

	#c("wgEncodeRegDnaseClusteredV2","wgEncodeRegTfbsClusteredV3","tfbsConsSites","gwasCatalog", "pubsBlat", "cosmic", "oreganno", "vistaEnhancers", "phastConsElements100way")

	tab.dt <- lapply(tables,getUCSCTable,genome=genome,cachedir=cachedir, sync=sync)
	names(tab.dt) <- tables

	togr <- function(x)
	{
		setnames(x,c("chrom","chromStart","chromEnd"),c("chr","start","end"))
		x$bin <- NULL
		x[,start:=start+1]
		makeGRanges(x)
	}
	tab.gr <- lapply(tab.dt,togr)
	names(tab.gr) <- tables
	return(tab.gr)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Generate feature sets based on CpG island, shore, and shelf regions
#'
#' Uses the "cpgIslandExt" table to generate shore (+/- 2kb from islands) and shelf (+/- 2kb from shores) regions. Will only function for genomes with this table available. Can be concatenated using c() with other tables, such as from getFeatures(), and provided as input for the "features" argument of goldmine().
#' @param genome The UCSC name specific to the genome of the query coordinates (e.g. "hg19", "hg18", "mm10", etc)
#' @param cachedir A path to a directory where a local cache of UCSC tables are stored. If equal to \code{NULL} (default), the data will be downloaded to temporary files and loaded on the fly. Caching is highly recommended to save time and bandwidth.
#' @export
getCpgFeatures <- function(genome, cachedir)
{
	island <- goldmine:::getFeatures("cpgIslandExt", genome, cachedir)[[1]]

	# Shores +/- 2kb from islands
	shore.r <- flank(island,2000,start=T,both=F)
	shore.l <- flank(island,2000,start=F,both=F)
	shore <- c(shore.r, shore.l)

	# Shelf +/- 2kb from shores
	shelf.r <- flank(shore.r,2000,start=T,both=F)
	shelf.l <- flank(shore.l,2000,start=F,both=F)
	shelf <- c(shelf.r, shelf.l)

	# Remove any conflicts
	# Island > Shore > Shelf
	shore <- setdiff(shore,island)
	shelf <- setdiff(shelf,island)
	shelf <- setdiff(shelf,shore)

	return(list(cpgIsland=island,cpgShore=shore,cpgShelf=shelf))
}
# --------------------------------------------------------------------

# ====================================================================



# ====================================================================
# Internal Functions

# --------------------------------------------------------------------
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
# --------------------------------------------------------------------

# --------------------------------------------------------------------
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
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# internal function to take two GRanges objects and output the % overlap of the query for each of the queries
# expands the output to include 0%s - should match row number of query.gr
calcPercentOverlap <- function(query.gr, subject.gr, sum.all=TRUE, report.bp=FALSE)
{
	# this function gives a percent without respect to strand (run twice on each strand separate for stranded results)
	# ignore strand if it was given in input GRanges
	strand(query.gr) <- "*"
	strand(subject.gr) <- "*"
	# need to set reduce the subject so overlapping ranges on different strand are now one range
	subject.gr <- reduce(subject.gr)

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
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# internal fucntion to take vectors of positions and a genome and output direct URLs to UCSC genome browser
getBrowserURLs <- function(input.gr, genome)
{
	# Give links to UCSC at the position of this region
	#http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr22%3A1-500
	paste("http://genome.ucsc.edu/cgi-bin/hgTracks?db=", genome, "&position=", seqnames(input.gr), "%3A", start(input.gr), "-", end(input.gr),sep="")
}
# --------------------------------------------------------------------
# ====================================================================

