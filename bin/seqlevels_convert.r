#!/usr/bin/env Rscript

library(GenomeInfoDb)
library(rtracklayer)
args <- commandArgs(trailingOnly=TRUE)

toUCSC = args[1]=="toUCSC"
inf = args[2]
## check file format
## if it is bigwig file
isBWF <- grepl("\\.(bw|bigwig)", inf, ignore.case=TRUE)
if(isBWF){## decrease the memory cost
    bwfile <- BigWigFile(inf)
    seqinfo <- seqinfo(bwfile)
    seqstyle <- seqlevelsStyle(seqinfo)
}else{
    data <- import(inf)
    seqstyle <- seqlevelsStyle(data)
}
readBWFile <- function(f, seqinfo){
    gr <- as(seqinfo, "GRanges")
    data <- GRanges()
    for(s in seq_along(gr)){
        dat <- import.bw(f, which = gr[s])
        dat <- coverage(dat, weight = dat$score)
        dat <- as(dat, "GRanges")
        dat <- dat[dat$score > 0] ## negative scores are not allowed
        data <- c(data, dat)
    }
    data <- coverage(data, weight = data$score)
    data <- as(data, "GRanges")
    data <- data[data$score > 0]
    return(data)
}
if(toUCSC){
    if(!"UCSC" %in% seqstyle){
        if(isBWF){
            data <- readBWFile(inf)
        }
        seqlevelsStyle(data) <- "UCSC"
        ## double check
        if(sum(grepl("^chr", seqlevels(data)))==0){
            ids <- grepl("^((\\d{1,2})|(IX|IV|V?I{0,3})|([XYMT]{1,2}))$", seqlevels(data))
            seqlevels(data)[ids] <- paste0("chr", seqlevels(data)[ids])
        }
        export(data, file.path(dirname(inf), paste0("UCSC.", basename(inf))))
    }else{
        file.copy(inf, file.path(dirname(inf), paste0("UCSC.", basename(inf))))
    }
}else{
    if(!"Ensembl" %in% seqstyle){
        if(isBWF){
            data <- readBWFile(inf)
        }
        seqlevelsStyle(data) <- "Ensembl"
        ## double check
        if(sum(grepl("^chr", seqlevels(data)))>0){
            ids <- grepl("^(chr)((\\d{1,2})|(IX|IV|V?I{0,3})|([XYMT]{1,2}))$", seqlevels(data))
            seqlevels(data)[ids] <- sub("chr", "", seqlevels(data)[ids])
        }
        export(data, file.path(dirname(inf), paste0("ENSEMBL.", basename(inf))))
    }else{
        file.copy(inf, file.path(dirname(inf), paste0("ENSEMBL.", basename(inf))))
    }
}

writeLines(as.character(packageVersion("GenomeInfoDb")), "GenomeInfoDb.version.txt")
writeLines(as.character(packageVersion("rtracklayer")), "rtracklayer.version.txt")
