#!/usr/bin/env Rscript

library(GenomeInfoDb)

args <- commandArgs(trailingOnly=TRUE)

inf = args[1] ## input file must be a bed file

data <- read.table(inf, nrows=1000, header=FALSE, quote=NULL, comment.char="#")

seqnames <- unique(as.character(data[, 1]))
seql <- "UCSC" %in% seqlevelsStyle(seqnames)

out <- "tmp.txt"

if(seql){
    writeLines("UCSC", out)
}else{
    writeLines("NOT_UCSC", out)
}

writeLines(as.character(packageVersion("GenomeInfoDb")), "GenomeInfoDb.version.txt")
