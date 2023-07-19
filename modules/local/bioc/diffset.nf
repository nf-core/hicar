process DIFFSET {
    tag "$bin_size"
    label 'process_low'

    conda "conda-forge::r-upsetr=1.0.3"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-upsetr:1.0.3--r3.3.2_1' :
        'biocontainers/r-upsetr:1.0.3--r3.3.2_1' }"

    input:
    tuple val(bin_size), path(peaks, stageAs: "peaks/*"), path(long_bedpe, stageAs: "long/*")
    val long_bedpe_postfix

    output:
    tuple val(bin_size), path("${prefix}/*")                                              , emit: diff
    tuple val(bin_size), val("$prefix"), path("${prefix}/setOperation.*") , optional: true, emit: anno
    path "versions.yml"                                                                   , emit: versions

    script:
    prefix   = task.ext.prefix ?: "setOperation_bin${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Nov. 11, 2022 call diffhic
    ## Copyright (c) 2022 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("UpSetR")
    versions <- c("NFCORE_HICAR:HICAR:DA:DIFFHIC:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
                    paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    OUTFOLDER = "$prefix"
    dir.create(OUTFOLDER, recursive = TRUE)

    ## get peaks
    pf <- dir("peaks", "bedpe", full.names = TRUE)
    header <- lapply(pf, read.delim, header=FALSE, nrow=1)
    peaks <- mapply(pf, header, FUN=function(f, h){
        hasHeader <- all(c("chr1", "start1", "end1",
                            "chr2", "start2", "end2") %in%
                            h[1, , drop=TRUE])
        .ele <- read.delim(f, header = hasHeader,
                            colClasses = rep("character", length(header)),
                            stringsAsFactors = FALSE)
        if(!hasHeader){
            colnames(.ele)[1:6] <- c("chr1", "start1", "end1",
                                "chr2", "start2", "end2")
        }
        .ele
    }, SIMPLIFY = FALSE)
    names(peaks) <- sub(".bedpe", "", gsub(".interactions.bedpe", "", basename(pf)))

    bedpe.ss <- lapply(peaks, function(.ele){
        with(.ele, paste(chr1, start1, end1, chr2, start2, end2, sep="\\t"))
    })
    symbols.all <- sort(unique(unlist(bedpe.ss)))

    symbols.mat <- matrix(0, nrow=length(symbols.all), ncol=length(bedpe.ss),
        dimnames=list(genes=symbols.all, groups=names(bedpe.ss)))
    for(i in seq_along(bedpe.ss)) symbols.mat[bedpe.ss[[i]], i] <- 1

    symbols.venn <- split(rownames(symbols.mat),
        apply(symbols.mat, 1, paste, collapse=""))
    names(symbols.venn) <-
        sapply(names(symbols.venn), function(.ele)
            paste(colnames(symbols.mat)[as.logical(as.numeric(strsplit(.ele, "")[[1]]))],
                collapse=".shared.with."))

    mapply(writeLines, symbols.venn,
        file.path(OUTFOLDER, paste0(names(symbols.venn), ".bedpe")))

    symbols.venn <- symbols.venn[!grepl("shared.with", names(symbols.venn))]
    if(length(symbols.venn)){
        mapply(writeLines, symbols.venn,
            file.path(OUTFOLDER,
                paste0("setOperation.", names(symbols.venn), ".only.bedpe")))
    }

    pdf("setOperation.upset.plot.pdf")
        upset(as.data.frame(symbols.mat))
    dev.off()
    """
}
