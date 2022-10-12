process MERGE_INTERACTIONS {
    tag "$fname"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bioconductor-rtracklayer=1.50.0" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'quay.io/biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(bin_size), path(fname)

    output:
    tuple val(bin_size), path("${bin_size}.bedpe")      , emit: interactions
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Oct. 12, 2022. Merge all the interaction files into one
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("GenomeInfoDb", "rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    inf = "$fname"
    ext = "${bin_size}.bedpe"
    data <- lapply(inf, import)
    data <- unlist(data)
    data <- reduce(data)
    export(data, ext)
    """
}
