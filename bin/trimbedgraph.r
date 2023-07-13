#!/usr/bin/env Rscript
#######################################################################
#######################################################################
## Trim bedgraph file by genome size
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
## This source code is licensed under the MIT license
## arg is the format of input files
#######################################################################
#######################################################################
args = commandArgs(TRUE)
parse_args <- function(options, args){
    lapply(options, function(.ele){
        n <- .ele[c(1, 2)]
        if(any(n %in% args)){
            if(.ele[3]=="logical"){
                list(arg=TRUE, ids=which(args %in% n))
            }else{
                id <- which(args %in% n)[1]
                x <- args[id+1]
                mode(x) <- .ele[3]
                list(arg=x, ids=c(id, id+1))
            }
        }else{
            NULL
        }
    })
}
option_list <- list("task_process"=c("--task_process", "-t", "character"),
                    "format"=c("--format", "-f", "character"),
                    "chrom_size"=c("--chrom_size", "-l", "character"))
opt <- parse_args(option_list, args)
to_be_removed <- unique(unlist(lapply(opt, function(.ele) .ele$ids)))
stopifnot(length(opt$format)>0)
format <- opt$format$arg
stopifnot(length(opt$chrom_size)>0)
genome <- read.delim(opt$chrom_size$arg, header=FALSE)
genome <- with(genome, GRanges(V1, IRanges(1, V2)))
fn <- args[-to_be_removed]
sl <- width(genome)
names(sl) <- as.character(seqnames(genome))
stopifnot(length(opt$task_process)>0)
pkgs <- c("rtracklayer")
versions <- paste0(opt$task_process$arg, ":")
for(pkg in pkgs){
    # load library
    library(pkg, character.only=TRUE)
    # parepare for versions.yml
    versions <- c(
        versions,
        paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
}
writeLines(versions, "versions.yml") # write versions.yml

# trim bed file
for(f in fn){
    d <- import(f, format=format, which=genome)
    if(tolower(format)=="wig"){
        d <- split(d, seqnames(d))
        d <- lapply(d, function(.d){
            w <- diff(end(.d))
            if(length(w)){
                w <- c(w, max(w, na.rm=TRUE))
                width(.d) <- w
            }
            .d
        })
        d <- unlist(GRangesList(d))
    }
    if(length(mcols(d)[, "score"])){
        d <- d[!is.na(mcols(d)[, "score"])]
    }
    seqlengths(d) <- sl[seqlevels(d)]
    d <- trim(d)
    seqlevels(d) <- sort(seqlevels(d), method='radix')
    d <- sort(d)
    n <- sub("\\..*?$", ".trimmed.bedgraph", f)
    export(d, n, format='bedGraph')
}
