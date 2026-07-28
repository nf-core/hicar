process MERGE_INTERACTIONS {
    tag "$fname"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(bin_size), path(fname)

    output:
    tuple val(bin_size), path("${prefix}.bedpe")        , emit: interactions
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "merged_${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Oct. 12, 2022. Merge all the interaction files into one
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

    inf = strsplit("$fname", "\\\\s+")[[1]]
    ext = "${prefix}.bedpe"
    data <- lapply(inf, readLines, n=1)
    inf <- inf[lengths(data)>0]
    if(length(inf)==0){
        data=GRanges()
    }else{
        data <- lapply(inf, import, format="BEDPE")
        data <- do.call(c, data)
        data <- unique(data)
        data <- sort(data)
        mcols(data)[, "name"] <- "*"
        mcols(data)[, "score"] <- 0
    }
    export(data, ext)
    """
}
