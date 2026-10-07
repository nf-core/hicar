process CALL_R1PEAK {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)
    path cut
    val pval

    output:
    tuple val(meta), path("*.narrowPeak")            , emit: peak
    path  "versions.yml"                             , emit: versions

    tuple val(meta), path("*.bed")                   , emit: bed
    tuple val(meta), path("*.bdg")                   , emit: bdg

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on DEC. 13, 2021 call R1 peaks
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    options(scipen=10)
    pkgs <- c("rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    ## options
    R1READS <- "$reads"
    CUT <- "$cut"
    r1reads <- read.delim(R1READS, header=FALSE)
    r1reads <- GRanges(r1reads[, 1], IRanges(r1reads[, 2], width=1), r1reads[, 6])
    digest <- import(CUT, format="BED")
    width(digest) <- 1
    ## peaks gap 5
    peaks <- reduce(digest, min.gapwidth=5L)
    mcols(peaks)[, "signalValue"] <- countOverlaps(peaks, r1reads, maxgap=5L)
    lambda <- mean(mcols(peaks)[, "signalValue"])
    p <- ppois(mcols(peaks)[, "signalValue"], lambda, lower.tail=FALSE)
    mcols(peaks)[, "pValue"] <- -10*log10(p)
    mcols(peaks)[, "qValue"] <- -10*log10(p.adjust(p, method="BH"))
    mcols(peaks)[, "score"] <- round((1-p)*100)
    peaks <- peaks[p<$pval]
    ## calculate the distance between digest sites
    dist <- distanceToNearest(digest)
    dist <- median(mcols(dist)[, "distance"])
    peaks <- promoters(peaks, upstream=dist, downstream=dist)
    peaks <- GenomicRanges::trim(peaks)
    rd <- reduce(peaks, min.gapwidth=dist, with.revmap=TRUE)
    revmap <- mcols(rd)[, "revmap"]
    l <- lengths(revmap)>1
    if(any(l)){
        rd_id <- unlist(revmap[l])
        rd_gp <- rep(seq_along(revmap[l]), lengths(revmap[l]))
        peaks_rd <- peaks[rd_id]
        peaks_rd_data <- split(as.data.frame(mcols(peaks_rd)), rd_gp)
        peaks_rd_data <- lapply(peaks_rd_data, FUN=colMeans)
        peaks_rd_data <- peaks_rd_data[order(as.numeric(names(peaks_rd_data)))]
        peaks_rd_data <- do.call(rbind, peaks_rd_data)
        peaks_rd <- rd[l]
        mcols(peaks_rd) <- peaks_rd_data
        peaks <- c(peaks[!l], peaks_rd)
        peaks <- sort(peaks)
    }
    if(any(start(peaks)<1)){
        start(peaks[start(peaks)<1]) <- 1
    }
    export(peaks, "${prefix}.bed")
    np <- paste(as.character(seqnames(peaks)), start(peaks)-1, end(peaks),
                ".", mcols(peaks)[, "score"],
                ".", mcols(peaks)[, "signalValue"],
                mcols(peaks)[, "pValue"],
                mcols(peaks)[, "qValue"],
                dist, sep="\t")
    writeLines(np, "${prefix}.narrowPeak")
    r1reads <- promoters(r1reads, upstream=75, downstream=75)
    cvg <- coverage(r1reads)
    export(cvg, "${prefix}_pileup.bdg", format="bedgraph")
    """
}
