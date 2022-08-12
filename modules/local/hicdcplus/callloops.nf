process HICDCPLUS_CALLLOOPS {
    tag "$bin_size"
    label 'process_high'
    label 'process_long'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::bioconductor-hicdcplus=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-hicdcplus:1.2.1--r41h619a076_0' :
        'quay.io/biocontainers/bioconductor-hicdcplus:1.2.1--r41h619a076_0' }"

    input:
    tuple val(meta), path(cool)
    path features

    output:
    tuple val(meta), path("*.interactions.txt")  , emit: interactions
    path "versions.yml"                          , emit: versions

    script:
    prefix   = task.ext.prefix ?: "diffhic_bin${bin_size}"
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on July 13, 2022 for HiC-DC+ to call signficant interactions
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################

    pkgs <- c("HiCDCPlus")
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
    """
}
