#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Aug. 24, 2021 enrichment analysis
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
c("optparse", "trackViewer")
# set libPath to pwd
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(trackViewer)
library(GenomicFeatures)
library(InteractionSet)
library(rtracklayer)
writeLines(as.character(packageVersion("trackViewer")), "trackViewer.version.txt")
writeLines(as.character(packageVersion("GenomicFeatures")), "GenomicFeatures.version.txt")
writeLines(as.character(packageVersion("InteractionSet")), "InteractionSet.version.txt")
writeLines(as.character(packageVersion("rtracklayer")), "rtracklayer.version.txt")

gap <- 1000
readwidth <- 150
maxEvent <- 25
if("optparse" %in% installed.packages()){
    library(optparse)
    option_list <- list(make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf", metavar="string"),
                        make_option(c("-s", "--chromsize"), type="character", default=NULL, help="filename of chrome size", metavar="string"),
                        make_option(c("-x", "--restrict"), type="character", default=NULL, help="filename of restrict cut", metavar="string"),
                        make_option(c("-r", "--resolution"), type="integer", default=NULL, help="resolution", metavar="integer"),
                        make_option(c("-d", "--gap"), type="integer", default=NULL, help="gap, default 2*resolution", metavar="integer"),
                        make_option(c("-l", "--readlength"), type="integer", default=NULL, help="reads length used to strengthen the signal, default resolution/10", metavar="integer"),
                        make_option(c("-m", "--maxevents"), type="integer", default=25, help="max events to plot", metavar="integer"),
                        make_option(c("-o", "--output"), type="character", default=".", help="output folder", metavar="string"),
                        make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"),
                        make_option(c("-e", "--events"), type="character", default=NULL, help="given events csv file, must be ginteractions file", metavar="string"))
    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
}else{
    args <- commandArgs(TRUE)
    parse_args <- function(options, args){
        out <- lapply(options, function(.ele){
            if(any(.ele[-3] %in% args)){
                if(.ele[3]=="logical"){
                    TRUE
                }else{
                    id <- which(args %in% .ele[-3])[1]
                    x <- args[id+1]
                    mode(x) <- .ele[3]
                    x
                }
            }
        })
    }
    option_list <- list("gtf"=c("--gtf", "-g", "character"),
                        "chromsize"=c("--chromsize", "-s", "character"),
                        "restrict"=c("--restrict", "-x", "character"),
                        "resolution"=c("--resolution", "-r", "integer"),
                        "gap"=c("--gap", "-d", "integer"),
                        "readlength"=c("--readlength", "-l", "integer"),
                        "maxevents"=c("--maxevents", "-m", "integer"),
                        "output"=c("--output", "-o", "character"),
                        "cores"=c("--cores", "-c", "integer"),
                        "events"=c("--events", "-e", "character"))
    opt <- parse_args(option_list, args)
}

gtf <- opt$gtf
resolution <- opt$resolution
gap <- 2 * resolution
readwidth <- ceiling(resolution/10)
output <- opt$output
chrom_size <- opt$chromsize
restrict_cut <- opt$restrict
if(!is.null(opt$gap)){
    gap <- as.numeric(opt$gap)
}
if(!is.null(opt$readlength)){
    readwidth <- as.numeric(opt$readlength)
}
if(!is.null(opt$maxevents)){
    maxEvent <- as.numeric(opt$maxevents)
}
cools <- dir(".", "mcool$", full.names = TRUE, recursive = TRUE)
pairfiles <- dir(".", ".unselected.pairs.gz$", full.names = TRUE, recursive = TRUE)
if(!is.null(opt$events)){
    evts <- read.csv(opt$events)
    if(!"fdr" %in% colnames(evts)){
        evts$fdr <- rep(1, nrow(evts))
    }
}else{
    evts <- dir(".", ".anno.csv$", full.names = TRUE, recursive = TRUE)
    if(length(evts)<1) stop("No events file.")
    evts <- lapply(evts, read.csv)
    evts <- do.call(rbind, evts)
}
stopifnot(nrow(evts)>0)
colnames(evts) <- tolower(colnames(evts))
stopifnot(all(c("chr1", "chr2", "start1", "start2", "end1", "end2", "fdr") %in%
    colnames(evts)))
evts <- evts[order(evts$fdr), , drop=FALSE]
if(nrow(evts)<1) stop("No events in files.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

## read file and summary counts
loadPairFire <- function(filenames){
    stopifnot(is.character(filenames))
    stopifnot(all(file.exists(filenames)))
    out <- list()
    for(fn in filenames){
        f <- gzfile(fn, open = "r")
        chunk <- readLines(f)
        close(f)
        chunk <- chunk[!grepl("^#", chunk)]
        if(length(chunk)){
            dt <- do.call(rbind, strsplit(chunk, "\\t"))
            dt <- unique(dt) # remove duplicates
            gi <- GInteractions(anchor1 = GRanges(dt[, 2], IRanges(as.numeric(dt[, 3]), width = readwidth), strand = dt[, 6]),
                                anchor2 = GRanges(dt[, 4], IRanges(as.numeric(dt[, 5]), width = readwidth), strand = dt[, 7]),
                                restrict1 = GRanges(dt[, 2], IRanges(as.numeric(dt[, 10]), as.numeric(dt[, 11]), names = dt[, 9])),
                                restrict2 = GRanges(dt[, 4], IRanges(as.numeric(dt[, 13]), as.numeric(dt[, 14]), names = dt[, 12])))
            out[[sub(".unselected.pairs.gz", "", basename(fn))]] <- gi
        }
    }
    return(out)
}
chunks <- loadPairFire(pairfiles)
readPairFile <- function(chunks, ranges){
    stopifnot(is(ranges, "GRanges"))
    out <- list()
    total <- list()
    for(fn in names(chunks)){
        chunk <- chunks[[fn]]
        if(length(chunk)){
            total[[fn]] <- length(chunk)
            out[[fn]] <-
                subsetByOverlaps(chunk, ranges = ranges, use.region="both")
        }
        rm(chunk)
        rm(gi)
    }
    ## split by group
    f <- sub("_REP\\d+", "", names(out))
    out <- split(out, f)
    giRbind <- function(a, b){
        GInteractions(anchor1 = c(first(a), first(b)),
                    anchor2 = c(second(a), second(b)),
                    restrict1 = c(a$restrict1, b$restrict1),
                    restrict2 = c(a$restrict2, b$restrict2))
    }
    out <- lapply(out, function(.ele){
        Reduce(giRbind, .ele)
    })
    total <- split(unlist(total), f)
    total <- sapply(total, sum)
    total <- median(total)/total
    return(list(gi=out, norm_factor=total))
}
restrict_pos <- import(restrict_cut, format = "BED")
seqlen <- read.delim(chrom_size, header=FALSE, row.names=1)
seql <- as.numeric(seqlen[, 1])
names(seql) <- rownames(seqlen)
txdb <- makeTxDbFromGFF(gtf)
gtf <- import(gtf, format = "GTF")
map <- gtf$gene_name
names(map) <- gtf$gene_id
grs <- with(evts, GInteractions(anchor1 = GRanges(chr1, IRanges(start1, end1)),
                                anchor2 = GRanges(chr2, IRanges(start2, end2)),
                                mode = "strict",
                                score = fdr))
grs$score <- 1-grs$score
grs <- unique(grs)
grs_narrow <- grs
regions(grs_narrow) <-
    resize(regions(grs), width = ceiling(width(regions(grs))/2), fix = "center")
link <- gi2track(grs_narrow)
setTrackStyleParam(link, "tracktype", "link")
setTrackStyleParam(link, "color", c("gray80", "yellow", "brown"))
reg <- reduce(GRanges(seqnames(first(grs)),
                        IRanges(start = start(first(grs)),
                                end = end(second(grs)))), min.gapwidth = gap)
gr1 <- unique(subsetByOverlaps(reg, grs))
gr1 <- gr1[seq.int(min(maxEvent, length(gr1)))]
prettyMax <- function(x){
    if(x<1) return(round(x, digits = 2))
    if(x<=5) return(ceiling(x))
    if(x<10) return(2*ceiling(x/2))
    if(x<100) return(10*ceiling(x/10))
    if(x<1000) return(50*ceiling(x/50))
    n <- 10^ceiling(log10(x))
    return(n*ceiling(x/n))
}
for(i in seq_along(gr1)){
    if(i < maxEvent){
        tryCatch({
            gr <- gr1[i]
            start(gr) <- start(gr) - 2* resolution
            end(gr) <- end(gr) + 2* resolution
            seqlevelsStyle(gr) <- seqlevelsStyle(txdb)[1]
            seqlengths(gr) <- seql
            gr <- GenomicRanges::trim(gr)
            chr_gr <- as(seqinfo(gr), "GRanges")[seqnames(gr)[1]]
            bait <- subsetByOverlaps(grs, gr, use.region="both")
            ## merge the bait by ends
            bait <- split(second(bait), as.character(first(bait)))
            bait <- lapply(bait, range)
            bait <- unlist(GRangesList(bait))
            bait <- GInteractions(anchor1 = parse2GRanges(names(bait)),
                                anchor2 = unname(bait))
            ## get the v4c
            info <- readPairFile(chunks, ranges = c(first(bait), second(bait)))
            total <- info$norm_factor
            names(cools) <- sub("\\d+\\.mcool$", "", basename(cools))
            gis <- lapply(cools, importGInteractions,
                        resolution = resolution,
                        format="cool",
                        ranges=gr, out = "GInteractions")
            gis <- mapply(gis, total[names(gis)], FUN=function(.ele, .total){
                .ele$score <- log2(.ele$score+1) * .total
                .ele
            })
            maxV <- max(sapply(gis, function(.ele) max(.ele$score, na.rm=TRUE)))
            maxV <- prettyMax(maxV)
            heat <- lapply(gis, function(.ele){
                .ele <- gi2track(.ele)
                setTrackStyleParam(.ele, "breaks", seq(from=0, to=maxV, by=maxV/10))
                setTrackStyleParam(.ele, "color", c("lightblue", "yellow", "red"))
                #setTrackStyleParam(.ele, "ylim", c(0, .5))
                .ele
            })
            names(heat) <- paste("valid pairs", names(heat))
            get_v4c <- function(bait1, bait2){
                v4c <- lapply(seq_along(bait1), function(.e){
                    vp1 <- mapply(info$gi, total[names(info$gi)], FUN=function(.ele, .total){
                        .ele <- subsetByOverlaps(.ele, bait1[.e], use.region="both")
                        .ele <- GRanges(coverage(c(first(.ele), second(.ele))))
                        .ele$score <- .ele$score*.total
                        new("track", dat=.ele,
                            type="data", format = "BED",
                            name = paste0(seqnames(bait1[.e]), ":",
                                        start(bait1[.e]), "-",
                                        end(bait1[.e])))
                    })
                    maxY <- sapply(vp1, FUN = function(.ele){
                        max(subsetByOverlaps(.ele$dat, bait2[.e])$score)
                    })
                    maxY <- prettyMax(max(maxY))
                    vp1 <- lapply(vp1, function(.ele){
                        setTrackStyleParam(.ele, "ylim", c(0, maxY))
                        setTrackStyleParam(.ele, "color", "gray80")
                        .ele
                    })
                })
                names(v4c) <- paste0(seqnames(bait1), ":",
                                    start(bait1), "-",
                                    end(bait1))
                v4c <- unlist(v4c, recursive = FALSE)
            }
            v4c_first <- get_v4c(first(bait), second(bait))
            v4c_second <- get_v4c(second(bait), first(bait))
            ids <- getGeneIDsFromTxDb(gr, txdb)
            genes <- geneTrack(ids, txdb, map[ids], asList=FALSE)
            tL <- trackList(genes, link, heat, v4c_first, v4c_second,
                            heightDist = c(1, 1, 2*length(heat),
                            2*length(v4c_first), 2*length(v4c_second)))
            names(tL)[2] <- "called links"
            pdf(file.path(output, paste0("event_", i, "_", seqnames(gr), ":", start(gr), "-", end(gr), ".pdf")),
                width = 9, height = length(tL))
            viewTracks(tL, gr=gr, autoOptimizeStyle = TRUE)
            dev.off()
        },
        error = function(e) message(e))
    }
}
