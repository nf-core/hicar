process SEQLEVELS_STYLE {
    tag "$bed"

    conda "bioconda::bioconductor-genomeinfodb=1.26.4"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0' :
        'biocontainers/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0' }"

    input:
    path bed

    output:
    stdout emit: seqlevels_style
    path "versions.yml"           , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created to detect the chromosome levels styles
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################

    library(GenomeInfoDb)
    versions <- c(
        "${task.process}:",
        paste("    GenomeInfoDb:", as.character(packageVersion("GenomeInfoDb"))))
    writeLines(versions, "versions.yml")

    inf = "$bed" ## input file must be a bed file

    data <- read.table(inf, nrows=1000, header=FALSE, quote=NULL, comment.char="#")

    seqnames <- unique(as.character(data[, 1]))
    seql <- "UCSC" %in% seqlevelsStyle(seqnames)

    if(seql){
        cat("UCSC")
    }else{
        cat("NOT_UCSC")
    }
    """
}
