#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on April. 29, 2021 call ChIPpeakAnno
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(ChIPpeakAnno)
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
writeLines(as.character(packageVersion("ChIPpeakAnno")), "ChIPpeakAnno.version.txt")

resList <- list()
args <- commandArgs(trailingOnly=TRUE)
gtf <- args[1]
pf <- file.path(args[2], "anno")
bin_size <- args[2]

detbl <- dir(".", "DEtable.*.csv|sig3Dinteractions.bedpe|peaks", recursive = TRUE, full.names = TRUE)
detbl <- detbl[!grepl("anno.csv", detbl)] ## in case of re-run

txdb <- makeTxDbFromGFF(gtf)
gtf <- import(gtf)
id2symbol <- function(gtf){
    if(is.null(gtf$gene_name)) return(NULL)
    x <- data.frame(id=gtf$gene_id, symbol=gtf$gene_name)
    x <- unique(x)
    x <- x[!duplicated(x$id), ]
    x <- x[!is.na(x$id), , drop=FALSE]
    if(nrow(x)==0) return(NULL)
    y <- x$symbol
    names(y) <- x$id
    y
}
id2symbol <- id2symbol(gtf)
anno <- toGRanges(txdb)
promoters <- promoters(anno, upstream=2000, downstream=500)
resList <- list()
peaks <- list()
promoterList <- list() # list to save the other sites of promoters

dir.create(pf, showWarnings = FALSE, recursive = TRUE)
for(det in detbl){
    if(grepl("csv$", det)) {
        DB <- read.csv(det)
    }else{
        if(grepl("peaks$", det)){
            DB <- read.table(det, header=TRUE)
        }else{
            DB <- read.delim(det)
        }
    }
    if(nrow(DB)<1) next
    rownames(DB) <- paste0("p", seq.int(nrow(DB)))
    DB.gr1 <- with(DB, GRanges(chr1, IRanges(start1, end1, name=rownames(DB))))
    DB.gr2 <- with(DB, GRanges(chr2, IRanges(start2, end2, name=rownames(DB))))
    # Annotation
    DB.anno1 <- annotatePeakInBatch(DB.gr1, AnnotationData = anno,
                                    output = "both",
                                    PeakLocForDistance = "middle",
                                    FeatureLocForDistance = "TSS",
                                    ignore.strand = TRUE)
    if(length(id2symbol)>0) DB.anno1$symbol[!is.na(DB.anno1$feature)] <- id2symbol[DB.anno1$feature[!is.na(DB.anno1$feature)]]
    DB.anno2 <- annotatePeakInBatch(DB.gr2, AnnotationData = anno,
                                    output = "both",
                                    PeakLocForDistance = "middle",
                                    FeatureLocForDistance = "TSS",
                                    ignore.strand = TRUE)
    if(length(id2symbol)>0) DB.anno2$symbol[!is.na(DB.anno2$feature)] <- id2symbol[DB.anno2$feature[!is.na(DB.anno2$feature)]]
    groupName <- sub(".sig3Dinteractions.bedpe|csv", "", basename(det))
    if(grepl("padj", det)){
        resList[[groupName]] <- c(DB.anno1, DB.anno2)
    }else{
        peaks[[groupName]] <- unique(c(DB.gr1, DB.gr2))
        ol1 <- findOverlaps(DB.gr1, promoters)
        ol2 <- findOverlaps(DB.gr2, promoters)
        promoterList[[groupName]] <- unique(c(DB.gr2[unique(queryHits(ol1))], DB.gr1[unique(queryHits(ol2))]))
    }
    # Summary the annotations
    DB.anno1 <- mcols(DB.anno1)
    DB.anno2 <- mcols(DB.anno2)
    DB.anno <- merge(DB.anno1, DB.anno2, by="peak",
                    suffixes = c(".anchor1",".anchor2"))
    DB <- cbind(DB[DB.anno$peak, ], DB.anno)
    pff <- file.path(pf, sub(".(csv|bedpe|peaks)", ".anno.csv", det))
    dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
    write.csv(DB, pff, row.names = FALSE)
}


if(packageVersion("ChIPpeakAnno")>="3.23.12"){
    if(length(resList)>0){
        if(is.list(resList)){
            resList <- GRangesList(resList[lengths(resList)>0])
        }
        out <- genomicElementDistribution(resList,
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                        # from 5' -> 3', fixed precedence 3' -> 5'
                                            breaks = c(-2000, -1000, -500, 0, 500),
                                            labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                    "upstream <500b", "TSS - 500b"),
                                            colors = c("#FFE5CC", "#FFCA99",
                                                    "#FFAD65", "#FF8E32")),
                                        plot = FALSE)

        ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".pdf")), plot=out$plot, width=9, height=9)
        ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".png")), plot=out$plot)
        out <- metagenePlot(resList, txdb)
        ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".pdf")), plot=out, width=9, height=9)
        ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".png")), plot=out)
    }
    if(length(peaks)>0){
        peaks <- GRangesList(peaks[lengths(peaks)>0])
        out <- genomicElementDistribution(peaks,
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                            # from 5' -> 3', fixed precedence 3' -> 5'
                                            breaks = c(-2000, -1000, -500, 0, 500),
                                            labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                    "upstream <500b", "TSS - 500b"),
                                            colors = c("#FFE5CC", "#FFCA99",
                                                    "#FFAD65", "#FF8E32")),
                                        plot = FALSE)

        ggsave(file.path(pf, paste0("genomicElementDistribuitonOfEachPeakList.", bin_size, ".pdf")), plot=out$plot, width=9, height=9)
        ggsave(file.path(pf, paste0("genomicElementDistribuitonOfEachPeakList.", bin_size, ".png")), plot=out$plot)

        out <- metagenePlot(peaks, txdb)
        ggsave(file.path(pf, paste0("metagenePlotToTSSOfEachPeakList.", bin_size, ".pdf")), plot=out, width=9, height=9)
        ggsave(file.path(pf, paste0("metagenePlotToTSSOfEachPeakList.", bin_size, ".png")), plot=out)

        if(length(peaks)<=5 && length(peaks)>1){
            ol <- findOverlapsOfPeaks(peaks)
            png(file.path(pf, paste0("vennDiagram.all.", bin_size, ".png")))
            vd <- makeVennDiagram(ol, connectedPeaks="keepAll")
            dev.off()
            write.csv(vd$vennCounts, file.path(pf, paste0("vennDiagram.all.", bin_size, ".csv")), row.names=FALSE)
        }
    }
    if(length(promoterList)>0){
        promoterList <- GRangesList(promoterList[lengths(promoterList)>0])
        out <- genomicElementDistribution(promoterList,
                                        TxDb = txdb,
                                        promoterRegion=c(upstream=2000, downstream=500),
                                        geneDownstream=c(upstream=0, downstream=2000),
                                        promoterLevel=list(
                                            # from 5' -> 3', fixed precedence 3' -> 5'
                                            breaks = c(-2000, -1000, -500, 0, 500),
                                            labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                    "upstream <500b", "TSS - 500b"),
                                            colors = c("#FFE5CC", "#FFCA99",
                                                    "#FFAD65", "#FF8E32")),
                                        plot = FALSE)

        ggsave(file.path(pf, paste0("genomicElementDistribuitonOfremoteInteractionPeaks.", bin_size, ".pdf")), plot=out$plot, width=9, height=9)
        ggsave(file.path(pf, paste0("genomicElementDistribuitonOfremoteInteractionPeaks.", bin_size, ".png")), plot=out$plot)

        if(length(promoterList)<=5 && length(promoterList)>1){
            ol <- findOverlapsOfPeaks(promoterList)
            png(file.path(pf, paste0("vennDiagram.remote.interaction.peak.with.promoters.all.", bin_size, ".png")))
            vd <- makeVennDiagram(ol, connectedPeaks="keepAll")
            dev.off()
            write.csv(vd$vennCounts, file.path(pf, paste0("vennDiagram.remote.interaction.peak.with.promoters.all.", bin_size, ".csv")), row.names=FALSE)
        }
    }
}
