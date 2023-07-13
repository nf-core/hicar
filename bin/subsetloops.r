#!/usr/bin/env Rscript
#######################################################################
#######################################################################
## Created on Oct. 26, 2022. Subset the interaction files by 1D peaks
## Copyright (c) 2022 Jianhong Ou (jianhong.ou@gmail.com)
## This source code is licensed under the MIT license
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
                    "prefix"=c("--prefix", "-a", "character"),
                    "peaks"=c("--peaks", "-p", "character"))
opt <- parse_args(option_list, args)
to_be_removed <- unique(unlist(lapply(opt, function(.ele) .ele$ids)))
stopifnot(length(opt$prefix)>0)
prefix <- opt$prefix$arg
stopifnot(length(opt$peaks)>0)
peaks <- read.delim(opt$peaks$arg, header=FALSE)
inf <- args[-to_be_removed]

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

ext <- paste0(prefix, inf)
peaks <- import(peaks)
data <- lapply(inf, import, format="BEDPE")
seqstyle1 <- seqlevelsStyle(peaks)
seqstyle2 <- seqlevelsStyle(first(data[[1]]))
if(length(intersect(seqstyle1, seqstyle2))==0){
    seqlevelsStyle(peaks) <- seqstyle2[1]
}
data <- lapply(data, function(.ele){
    keep <- countOverlaps(first(.ele), subject=peaks, ignore.strand=TRUE) > 0 |
        countOverlaps(second(.ele), subject=peaks, ignore.strand=TRUE)>0
    .ele[keep]
})

dir.create("sub_loop")
mapply(export, data, file.path("sub_loop", ext), format="BEDPE")

data2 <- lapply(data, function(.ele){
    keep <- countOverlaps(query=peaks, subject=first(.ele), ignore.strand=TRUE) > 0 |
        countOverlaps(query=peaks, subject=second(.ele), ignore.strand=TRUE) > 0
    peaks[keep]
})

dir.create("sub_peak")
mapply(export, data2,
    file.path("sub_peak", paste0(ext, ".bed")),
    format="BED")
