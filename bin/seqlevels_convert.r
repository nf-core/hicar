#!/usr/bin/env Rscript

library(GenomeInfoDb)
library(rtracklayer)
args <- commandArgs(trailingOnly=TRUE)

toUCSC = args[1]=="toUCSC"
inf = args[2]

data <- import(inf)

if(toUCSC){
    if(!"UCSC" %in% seqlevelsStyle(data)){
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
    if(!"Ensembl" %in% seqlevelsStyle(data)){
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
