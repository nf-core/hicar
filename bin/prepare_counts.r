#!/usr/bin/env Rscript
#######################################################################
#######################################################################
## Created on Aug. 24, 2021 count reads for peak filtering
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################

pkgs <- c("rtracklayer", "InteractionSet", "Biostrings", "Rsamtools")
for(pkg in pkgs){
    library(pkg, character.only=TRUE)
    writeLines(as.character(packageVersion(pkg)), paste0(pkg, ".version.txt"))
}

OUTPUT <- "counts.csv"

if("optparse" %in% installed.packages()){
    library(optparse)
    option_list <- list(make_option(c("-a", "--r1peak"), type="character", default=NULL, help="filename of r1 peak", metavar="string"),
                        make_option(c("-b", "--r2peak"), type="character", default=NULL, help="filename of r2 peak", metavar="string"),
                        make_option(c("-x", "--restrict"), type="character", default=NULL, help="filename of restrict cut", metavar="string"),
                        make_option(c("-p", "--pairs"), type="character", default=NULL, help="folder of valid distal pairs", metavar="string"),
                        make_option(c("-m", "--mappability"), type="character", default=NULL, help="mappability file", metavar="string"),
                        make_option(c("-o", "--output"), type="character", default="counts.csv", help="output folder", metavar="string"),
                        make_option(c("-f", "--fasta"), type="character", default=NULL, help="genome fasta file", metavar="string"))
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
    option_list <- list("r1peak"=c("--r1peak", "-a", "character"),
                        "r2peak"=c("--r2peak", "-b", "character"),
                        "restrict"=c("--restrict", "-x", "character"),
                        "mappability"=c("--mappability", "-m", "character"),
                        "pairs"=c("--pairs", "-p", "character"),
                        "output"=c("--output", "-o", "character"),
                        "fasta"=c("--fasta", "-f", "character"))
    opt <- parse_args(option_list, args)
}
for(i in c("r1peak", "r2peak", "restrict", "mappability", "pairs", "fasta")){
    if(is.null(opt[[i]])){
        stop(i, " is required!")
    }
}
if(!is.null(opt$output)){
    OUTPUT <- opt$output
}
FASTA <- opt$fasta
CUT <- opt$restrict
MAPPABILITY <- opt$mappability
pairs <- dir(opt$pairs, full.names=TRUE)
R1PEAK <- import(opt$r1peak)
R2PEAK <- import(opt$r2peak)
## split by chromsome
R1PEAK <- split(R1PEAK, seqnames(R1PEAK))
R2PEAK <- split(R2PEAK, seqnames(R2PEAK))
R1PEAK <- R1PEAK[lengths(R1PEAK)>0]
R2PEAK <- R2PEAK[lengths(R2PEAK)>0]
chromosomes <- intersect(names(R1PEAK), names(R2PEAK))
if(length(chromosomes)==0){
    stop("no valid data in same chromosome.")
}
## loading data


countByOverlaps <- function(gi, pairs){
    counts_tab <- 0
    counts_r2 <- 0
    for(ps in pairs){
        ps <- read.delim(ps, comment.char="#", header=FALSE)
        ps <- with(ps, GInteractions(GRanges(V2, IRanges(V3, width=100)),
                                    GRanges(V4, IRanges(V5, width=100))))
        counts_tab <- counts_tab + countOverlaps(gi, ps, use.region="same")
        counts_r2 <- counts_r2 + countOverlaps(second(gi), second(ps))
    }
    list(count=counts_tab, shortCount=counts_r2)
}

getMscore <- function(mscore, gr){
    vw <- Views(mscore, unique(gr))
    sc <- viewMeans(vw)
    vs <- ranges(vw)
    vs <- as(vs, "GRanges")
    vs$score <- unlist(sc, use.names=FALSE)
    ol <- findOverlaps(gr, vs, type="equal")
    ol <- ol[!duplicated(queryHits(ol))]
    score <- vs[subjectHits(ol)]$score[match(seq_along(gr), queryHits(ol))]
    score[is.na(score)] <- 0
    score
}

### load counts
gis <- NULL
for(chrom in chromosomes){
    r1peak <- R1PEAK[[chrom]]
    r2peak <- R2PEAK[[chrom]]
    peak_pair <- expand.grid(seq_along(r1peak), seq_along(r2peak))
    gi <- GInteractions(r1peak[peak_pair[, 1]], r2peak[peak_pair[, 2]])
    cnt <- countByOverlaps(gi, pairs)
    gi$count <- cnt$count
    gi$shortCount <- cnt$shortCount
    gi <- gi[gi$count>0]
    if(length(gi)>0){
        gis <- c(gis, gi)
    }
}
gis <- do.call(c, gis)
### load gc content
fa <- FaFile(file=FASTA)
gc1 <- letterFrequency(getSeq(fa, first(gis)), letters="CG", as.prob=TRUE)
gc2 <- letterFrequency(getSeq(fa, second(gis)), letters="CG", as.prob=TRUE)
gis$gc <- gc1 * gc2 + 1e-9
### load enzyme cut number in R1
cut <- import(CUT)
start(cut) <- end(cut)
gis$cut <- countOverlaps(first(gis), cut)+0.1
### load mapping score
mscore <- import(MAPPABILITY, as="RleList")
m1 <- getMscore(mscore, first(gis))
m2 <- getMscore(mscore, second(gis))
gis$mappability <- m1*m2 + 1e-6
### get distance of the anchors
gis$dist <- distance(first(gis), second(gis))+1
gis$length <- width(first(gis))*width(second(gis))

mm <- log(as.data.frame(mcols(gis)[, c("length", "cut", "gc", "mappability", "dist", "shortCount")]))
mm <- cbind(as.data.frame(first(gis)), as.data.frame(second(gis)), gis$count, mm)
colnames(mm) <- c("chr1", "start1", "end1", "width1", "strand1",
                "chr2", "start2", "end2", 'width2', "strand2",
                "count", "logl", "logn", "loggc", "logm", "logdist", 'logShortCount')
mm$strand1 <- NULL
mm$strand2 <- NULL
write.csv(mm, OUTPUT, row.names=FALSE)
