#!/usr/bin/env Rscript
#######################################################################
#######################################################################
## Created on Oct. 2021 convert pairs.gz to bam file for visualization
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################

library(Rsamtools)
writeLines(as.character(packageVersion("Rsamtools")), "Rsamtools.version.txt")
library(InteractionSet)
writeLines(as.character(packageVersion("InteractionSet")), "InteractionSet.version.txt")

args <- commandArgs(TRUE)
peaks <- args[1]
pairs <- dir(".", "unselected.pairs.gz")
cons<- sub("unselected.pairs.gz", "bam", pairs)

scanPairHeader <- function(file){
    con <- gzcon(file(file, open="rb"))
    on.exit(close(con))
    header <- NULL
    while(length(line <- readLines(con, n = 1))){
        if(grepl("^#", line)){
            header <- c(header, line)
        }else{
            return(header)
        }
    }
}

dataInPeak <- function(data, peak){
    dgi <- with(data, GInteractions(GRanges(V2, IRanges(V3, width=150)),
                                    GRanges(V4, IRanges(V5, width=150))))
    data[countOverlaps(dgi, peak, use.region="same")>0, , drop=FALSE]
}

exportBamFile <- function(header, data, con){
    sam_path <- sub("bam$", "sam", con, ignore.case = TRUE)
    if(sam_path==con){
        sam_path <- paste0(con, ".sam")
    }
    sam_con <- file(sam_path, "w")
    on.exit(close(sam_con))
    data <- data[data[, 2]==data[, 4], , drop=FALSE]
    #aln1 <- data[, c(1, 6, 2, 3)]
    #aln2 <- data[, c(1, 7, 4, 5)]
    #colnames(aln1) <- colnames(aln2) <- c("name", "strand", "chr", "pos")
    #aln1$name <- paste0(aln1$name, "/1")
    #aln2$name <- paste0(aln2$name, "/2")
    #getFlag <- function(s1, s2, first=TRUE){
    #    3+ifelse(s1=="-", 16, 0) + ifelse(s2=="-", 32, 0) + ifelse(first, 64, 128)
    #}
    #aln1 <- cbind(aln1, mrnm=aln2$chr, mpos=aln2$pos, flag=getFlag(aln1$strand, aln2$strand, TRUE))
    #aln2 <- cbind(aln2, mrnm=aln1$chr, mpos=aln1$pos, flag=getFlag(aln2$strand, aln1$strand, FALSE))
    #data <- rbind(aln1, aln2)
    #data <- data[order(rep(seq.int(nrow(data)/2), 2)), ]
    #rm(aln1, aln2)
    #aln <- paste(data$name,
    #            data$flag,
    #            data$chr,
    #            data$pos,
    #            "50",
    #            "100M",
    #            data$mrnm,
    #            data$mpos,
    #            abs(data$pos - data$mpos),
    #            "*",
    #            "*",
    #            sep = "\t")
    colnames(data) <- c("name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")
    data$flag <- ifelse(data$strand1=="-", 16, 0)
    data$pos <- rowMins(as.matrix(data[, c("pos1", "pos2")]))
    data$isize <- abs(data$pos1 - data$pos2)
    data$cigar <- paste0("100M", data$isize, "N100M")
    aln <- paste(data$name, data$flag, data$chr1, data$pos, "50", data$cigar, "*", 0, data$isize+201, "*", "*", sep = "\t")

    writeLines(header, sam_con)
    writeLines(aln, sam_con)
    close(sam_con)
    on.exit()
    si <- do.call(rbind, strsplit(header, "\\t"))
    si <- as.numeric(sub("LN:", "", si[, 3]))
    si <- si[!is.na(si)]
    if(length(si)){
        si <- any(si>536870912)
    }else{
        si <- TRUE
    }
    if(si){
        bam <- asBam(sam_path, sub(".bam$", "", con, ignore.case = TRUE),
                    overwrite = TRUE, indexDestination = FALSE)
    }else{
        bam <- asBam(sam_path, sub(".bam$", "", con, ignore.case = TRUE),
                    overwrite = TRUE, indexDestination = TRUE)
    }
    unlink(sam_path)
    invisible(bam)
}

peaks <- read.csv(peaks)
peaks <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)),
                                    GRanges(chr2, IRanges(start2, end2))))

null <- mapply(function(file, con){
    header <- scanPairHeader(file)
    header <- header[grepl("#samheader: @SQ", header)]
    header <- sub("#samheader: ", "", header)

    data <- read.table(file,
                    colClasses=c("character", "character",
                                "integer", "character",
                                "integer", "character",
                                "character",
                                rep("NULL", 7)))
    data <- data[!duplicated(data[, -1]), , drop=FALSE]
    data <- dataInPeak(data, peaks)
    exportBamFile(header, data, con)
}, pairs, cons)
