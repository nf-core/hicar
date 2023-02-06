process BEDGRAPH_TRIM {
    tag "$bedgraph"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'quay.io/biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.trimmed.bedgraph")         , emit: bedgraph
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: 'bedGraph'
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Trim bedgraph file by genome size
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    ## arg is the format of input files
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

    format="$args"
    fn <- strsplit("$bedgraph", " ")[[1]]
    genome <- read.delim("$sizes", header=FALSE)
    genome <- with(genome, GRanges(V1, IRanges(1, V2)))
    for(f in fn){
        d <- import(f, format=format, which=genome)
        if(tolower(format)=="wig"){
            w <- diff(end(d))
            if(length(w)){
                w <- c(w, max(w, na.rm=TRUE))
                width(d) <- w
            }
        }
        if(length(mcols(d)[, "score"])){
            d <- d[!is.na(mcols(d)[, "score"])]
        }
        sl <- width(genome)
        names(sl) <- as.character(seqnames(genome))
        seqlengths(d) <- sl
        d <- trim(d)
        n <- sub("\\\\..*?\$", ".trimmed.bedgraph", f)
        export(d, n, format='bedGraph')
    }
    """
}
