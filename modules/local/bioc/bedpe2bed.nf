process BEDPE2BED {
    tag "$fname"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(meta), path(fname)

    output:
    tuple val(meta), path("${prefix}_regions.bed")      , emit: bed
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "$meta.id"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created to convert bedpe format to bed format
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    pf <- dir(".", "bedpe", full.names = TRUE)
    header <- lapply(pf, read.delim, header=FALSE, nrow=1)
    peaks <- mapply(pf, header, FUN=function(f, h){
        hasHeader <- all(c("chr1", "start1", "end1", "chr2", "start2", "end2") %in%
                        h[1, , drop=TRUE])
        .ele <- read.delim(f, header = hasHeader)
        if(!hasHeader){
            colnames(.ele)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
            ## bedpe, [start, end)
            .ele\$start1 <- .ele\$start1+1
            .ele\$start2 <- .ele\$start2+1
        }
        .ele
    }, SIMPLIFY = FALSE)
    ### reduce the peaks
    peaks <- unique(do.call(rbind, peaks)[, c("chr1", "start1", "end1",
        "chr2", "start2", "end2")])
    peaks <- with(peaks, c(GRanges(chr1, IRanges(start1, end1)),
        GRanges(chr2, IRanges(start2, end2))))
    peaks <- sort(unique(peaks))
    export(peaks, "${prefix}_regions.bed")
    """
}
