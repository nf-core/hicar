process SEQLEVELS_STYLE {
    tag "$bed"

    conda (params.enable_conda ? "bioconda::bioconductor-genomeinfodb=1.26.4" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0' }"

    input:
    path bed

    output:
    stdout emit: seqlevels_style
    path "versions.yml"           , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

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
