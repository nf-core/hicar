#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Sept. 10, 2021 stats for ATAC reads
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

args <- commandArgs(trailingOnly=TRUE)
gtf <- args[1]

readsFiles <- dir(".", "merged.ATAC.bed.gz") ## postfix from mergedreads.nf
peaksFiles <- dir(".", "narrowPeak|broadPeak") ## output from macs2

names(readsFiles) <- sub(".merged.ATAC.bed.gz", "", readsFiles)
names(peaksFiles) <- sub("_peaks.*Peak", "", peaksFiles)

N <- intersect(names(readsFiles), names(peaksFiles))
stopifnot("no peak and signal pairs"=length(N)>0)

peaksFiles <- peaksFiles[N]
readsFiles <- readsFiles[N]

library(rtracklayer)
library(GenomicFeatures)
library(ATACseqQC)
writeLines(as.character(packageVersion("ATACseqQC")), "ATACseqQC.version.txt")
writeLines(as.character(packageVersion("GenomicFeatures")), "GenomicFeatures.version.txt")
writeLines(as.character(packageVersion("rtracklayer")), "rtracklayer.version.txt")

## import reads
readls <- lapply(readsFiles, function(f){
    reads <- NULL
    con <- file(f, "r")
    on.exit(close(con))
    while(length(line <- readLines(con, n = 1e6))){
        line <- sub("\\t\\*\\t\\*", "", line)
        tab <- table(line)
        tab <- as.data.frame(tab)
        reads <- rbind(reads, tab)
        reads <- aggregate(Freq ~ line, data = reads, FUN = sum, drop = FALSE)
    }
    reads <- cbind(reads, do.call(rbind, strsplit(as.character(reads[, 1]), "\\t")))
    reads <- GRanges(reads[, 3],
                    IRanges(as.numeric(reads[, 4])+1, as.numeric(reads[, 4])+150),
                    strand = reads[, 6],
                    score = reads[, 2])
})

## import peaks
peakls <- lapply(peaksFiles, import)

txdb <- makeTxDbFromGFF(gtf)
txs <- exons(txdb)

stats <- mapply(peakls, readls, FUN = function(peaks, reads){
    ## calculate FRiP score (fraction of reads in peaks), must over 1%
    readsInPeaks <- subsetByOverlaps(reads, peaks)
    FRiP <- 100*sum(readsInPeaks$score)/sum(reads$score)

    gal <- as(reads, "GAlignments")
    gal <- rep(gal, mcols(gal)$score)
    ## calculate Transcription Start Site Enrichment Score
    tsse <- TSSEscore(gal, txs)

    ## Promoter/Transcript body (PT) score
    pt <- PTscore(gal, txs)
    pt <- pt[(pt$promoter>0 | pt$transcriptBody>0) & pt$PT_score!=0]
    promoterEnriched <- table(pt$PT_score>0)
    names(promoterEnriched) <-
        c("FALSE"="bodyEnrich", "TRUE"="promoterEnrich")[names(promoterEnriched)]
    promoterEnriched <-
        c(promoterEnriched,
            prop.test=prop.test(cbind(table(pt$PT_score>0), c(50, 50)))$p.value)
    list(tsse_FRiP = c(TSSEscore=tsse$TSSEscore, FRiP=FRiP, promoterEnriched),
        tsseValues = tsse$values)
}, SIMPLIFY = FALSE)

tsse_FRiP <- do.call(rbind, lapply(stats, function(.ele) .ele$tsse_FRiP))
write.csv(tsse_FRiP, "TSSEscore_FRiP.csv")

tsse <- do.call(rbind, lapply(stats, function(.ele) .ele$tsseValues))
colnames(tsse) <- 100*(-9:10-.5)
write.csv(tsse, "aggregateTSSEscoreToTSS.csv")
