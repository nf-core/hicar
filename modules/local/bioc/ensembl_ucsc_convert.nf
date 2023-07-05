process ENSEMBL_UCSC_CONVERT {
    tag "$fname"
    label 'process_medium'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(bin_size), path(fname)

    output:
    tuple val(bin_size), path("{UCSC,ensembl}.${fname}"), emit: tab
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created to convert the chromosome levels styles
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
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

    toUCSC = "$args"=="toUCSC"
    inf = "$fname"
    ## check file format
    ## if it is bigwig file
    isBWF <- grepl("\\\\.(bw|bigwig)", inf, ignore.case=TRUE)
    if(isBWF){## decrease the memory cost
        bwfile <- BigWigFile(inf)
        seqinfo <- seqinfo(bwfile)
        seqstyle <- seqlevelsStyle(seqinfo)
    }else{
        data <- import(inf)
        seqstyle <- seqlevelsStyle(data)
    }
    readBWFile <- function(f, seqinfo){
        gr <- as(seqinfo, "GRanges")
        data <- GRanges()
        for(s in seq_along(gr)){
            dat <- import.bw(f, which = gr[s])
            dat <- coverage(dat, weight = dat\$score)
            dat <- as(dat, "GRanges")
            dat <- dat[dat\$score > 0] ## negative scores are not allowed
            data <- c(data, dat)
        }
        data <- coverage(data, weight = data\$score)
        data <- as(data, "GRanges")
        data <- data[data\$score > 0]
        return(data)
    }
    if(toUCSC){
        if(!"UCSC" %in% seqstyle){ ## convert to UCSC style
            if(isBWF){
                data <- readBWFile(inf, seqinfo)
            }
            seqlevelsStyle(data) <- "UCSC"
            ## double check
            if(sum(grepl("^chr", seqlevels(data)))==0){
                ids <- grepl("^((\\\\d{1,2})|(IX|IV|V?I{0,3})|([XYMT]{1,2}))\$", seqlevels(data))
                seqlevels(data)[ids] <- paste0("chr", seqlevels(data)[ids])
            }
            export(data, file.path(dirname(inf), paste0("UCSC.", basename(inf))))
        }else{
            file.copy(inf, file.path(dirname(inf), paste0("UCSC.", basename(inf))))
        }
    }else{
        if(!"Ensembl" %in% seqstyle){## convert to Ensembl style
            if(isBWF){
                data <- readBWFile(inf, seqinfo)
            }
            seqlevelsStyle(data) <- "Ensembl"
            ## double check
            if(sum(grepl("^chr", seqlevels(data)))>0){
                ids <- grepl("^(chr)((\\\\d{1,2})|(IX|IV|V?I{0,3})|([XYMT]{1,2}))\$", seqlevels(data))
                seqlevels(data)[ids] <- sub("chr", "", seqlevels(data)[ids])
            }
            export(data, file.path(dirname(inf), paste0("ENSEMBL.", basename(inf))))
        }else{
            file.copy(inf, file.path(dirname(inf), paste0("ENSEMBL.", basename(inf))))
        }
    }
    """
}
