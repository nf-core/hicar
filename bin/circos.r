#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Oct. 16, 2021 prepare data for circos
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################

interaction <- ""
chromsize <- ""
gtf <- ""
ucscname <- ""
outfolder <- "circos"

pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

options(scipen=10)
library(rtracklayer)
writeLines(as.character(packageVersion("rtracklayer")), "rtracklayer.version.txt")

if("optparse" %in% installed.packages()){
    library(optparse)
    option_list <- list(make_option(c("-i", "--interaction"), type="character", default=NULL, help="interaction bedpe file", metavar="string"),
                        make_option(c("-g", "--gtf"), type="character", default=NULL, help="annotation gtf file", metavar="string"),
                        make_option(c("-c", "--chromsize"), type="character", default=NULL, help="filename of chromosome size", metavar="string"),
                        make_option(c("-u", "--ucscname"), type="character", default=NULL, help="ucsc annotation name", metavar="string"))
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
    option_list <- list("interaction"=c("--interaction", "-i", "character"),
                        "gtf"=c("--gtf", "-g", "character"),
                        "chromsize"=c("--chromsize", "-c", "character"),
                        "ucscname"=c("--ucscname", "-u", "character"))
    opt <- parse_args(option_list, args)
}
for(i in c("interaction", "gtf", "chromsize", "ucscname")){
    if(is.null(opt[[i]])){
        stop(i, " is required")
    }else{
        assign(i, opt[[i]])
    }
}

dir.create(outfolder, showWarnings = FALSE)

pe <- import(interaction, format="BEDPE")
seqlevelsStyle(first(pe)) <- seqlevelsStyle(second(pe)) <- "UCSC"
pes <- pe[order(mcols(pe)$score, decreasing=TRUE)]
pes_cis <- pes[seqnames(first(pe))==seqnames(second(pe))]
pes_trans <- pes[seqnames(first(pe))!=seqnames(second(pe))]
pes <- sort(c(pes_cis[seq.int(min(1e4, length(pes_cis)))],
            pes_trans[seq.int(min(1e4, length(pes_trans)))])) ## keep top 10K links only. otherwise hard to plot.
out <- as.data.frame(pes)
scores <- sqrt(range(mcols(pe)$score)/10)
scores <- c(floor(scores[1]), ceiling(scores[2]))
cid <- cut(sqrt(mcols(pes)$score/10), breaks = seq(scores[1], scores[2]))
levels(cid) <- seq_along(levels(cid))
out <- cbind(out[, c("first.seqnames", "first.start", "first.end",
                    "second.seqnames", "second.start", "second.end")],
             thickness=paste0("thickness=", cid))

write.table(out, file.path(outfolder, "link.txt"),
    quote=FALSE, col.names=FALSE, row.names=FALSE,
    sep=" ")


chromsize <- read.delim(chromsize, header=FALSE)
seqlevelsStyle(chromsize[, 1]) <- "UCSC"
chromsize <- cbind("chr", "-", chromsize[, c(1, 1)], 0, chromsize[, 2],
                    paste0("chr", sub("^chr", "", chromsize[, 1])))
colnames(chromsize) <- c("chr", "-", "chrname", "chrlabel",
                        0, "chrlen", "chrcolor")
chromsize <- chromsize[order(as.numeric(sub("chr", "", chromsize[, "chrname"])),
                        chromsize[, "chrname"]), , drop=FALSE]
write.table(chromsize, file.path(outfolder, "karyotype.tab"),
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")

getScore <- function(seql, rg){
    gtile <- tileGenome(seqlengths = seql, tilewidth = 1e5)
    gtile <- unlist(gtile)
    gtile <- split(gtile, seqnames(gtile))
    sharedChr <- intersect(names(rg), names(gtile))
    rg <- rg[sharedChr]
    gtile <- gtile[names(rg)]
    vw <- Views(rg, gtile)
    vm <- viewSums(vw, na.rm = TRUE)
    stopifnot(identical(names(vm), names(gtile)))
    gtile <- unlist(gtile, use.names = FALSE)
    gtile$score <- unlist(vm, use.names = FALSE)
    gtile
}
rg <- coverage(c(first(pe), second(pe)))
seql <- chromsize$chrlen
names(seql) <- chromsize$chrname
gtile <- getScore(seql, rg)
rg <- as.data.frame(gtile)
rg <- rg[rg$score>0, c("seqnames", "start", "end", "score"), drop=FALSE]
write.table(rg, file.path(outfolder, "hist.link.txt"),
            quote=FALSE, col.names=FALSE, row.names=FALSE,
            sep=" ")
seqn <- sort(as.character(unique(rg$seqnames)))[1]
labelA <- labelB <- c(seqn, 0, 5000)
labelA <- c(labelA, "interaction-density")
labelB <- c(labelB, "exon-density")
writeLines(labelA, file.path(outfolder, "labelA.txt"), sep=" ")
writeLines(labelB, file.path(outfolder, "labelB.txt"), sep=" ")

gtf <- import(gtf)
gtf <- gtf[gtf$type %in% "exon"]
seqlevelsStyle(gtf) <- "UCSC"
rg <- coverage(gtf)
gtile <- getScore(seql, rg)
rg <- as.data.frame(gtile)
rg <- rg[rg$score>0, c("seqnames", "start", "end", "score"), drop=FALSE]
write.table(rg, file.path(outfolder, "hist.exon.txt"),
            quote=FALSE, col.names=FALSE, row.names=FALSE,
            sep=" ")
tryCatch({
    session <- browserSession()
    genome(session) <- ucscname
    ideo <- getTable(ucscTableQuery(session, table="cytoBandIdeo"))
    ideo <- ideo[ideo$chrom %in% chromsize$chrname, , drop=FALSE]
    ideo <- data.frame("band", ideo$chrom, ideo$name, ideo$name,
                ideo$chromStart, ideo$chromEnd, ideo$gieStain)
    colnames(ideo) <- colnames(chromsize)
    ideo <- rbind(chromsize, ideo)
    for(i in c("chrname", "chrlabel")){
        ideo[ideo[, i]=="", i] <- make.names(ideo[ideo[, i]=="", 2], unique=TRUE)
    }
    write.table(ideo, file.path(outfolder, "karyotype.tab"),
                quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
}, error=function(.e) message(.e))
