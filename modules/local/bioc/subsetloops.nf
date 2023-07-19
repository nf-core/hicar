process BIOC_SUBSETLOOPS {
    tag "$mergedloops"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    path peaks
    tuple val(bin_size), path(mergedloops)

    output:
    tuple val(bin_size), path("sub_loop/${prefix}*")    , emit: loops
    tuple val(bin_size), path("sub_peak/${prefix}*")    , emit: bed
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'subset_'
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Oct. 26, 2022. Subset the interaction files by 1D peaks
    ## Copyright (c) 2022 Jianhong Ou (jianhong.ou@gmail.com)
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

    inf <- strsplit("$mergedloops", "\\\\s+")[[1]]
    ext <- paste0("$prefix", inf)
    peaks <- import("$peaks")
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
    """
}
