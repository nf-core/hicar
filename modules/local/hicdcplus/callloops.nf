process HICDCPLUS_CALLLOOPS {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::bioconductor-hicdcplus=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-hicdcplus:1.2.1--r41h619a076_0' :
        'quay.io/biocontainers/bioconductor-hicdcplus:1.2.1--r41h619a076_0' }"

    input:
    tuple val(meta), path(cool), val(bin_size), path(features)
    path chromsize

    output:
    tuple val(meta), path("*.interactions.txt")  , emit: interactions
    path "versions.yml"                          , emit: versions

    script:
    prefix   = task.ext.prefix ?: "$meta.id"
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on July 13, 2022 for HiC-DC+ to call signficant interactions
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################

    pkgs <- c("HiCDCPlus", "dplyr")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    ## step 1, generate features
    bintolen <- data.table::fread("$features")
    bintolen <- bintolen %>%
        tidyr::separate(.data$bins, c("chr", "start", "end"), "-") %>%
            dplyr::mutate(start = floor(as.numeric(.data$start)/1000) *
                1000, end = as.numeric(.data$end))
    binsize <- $bin_size
    seqlen <- read.delim("$chromsize")

    gi_list <- generate_bintolen_gi_list(bintolen_path="$features")

    """
}
