process PAIR2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(peak), path(pairs)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai")    , emit: bam
    path "versions.yml"                                  , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Oct. 2021 convert pairs.gz to bam file for visualization
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("Rsamtools", "InteractionSet", "rhdf5")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    peaks <- "$peak"
    pairs <- dir(".", "h5\$")

    ## load header
    getHeader <- function(file){
        header <- h5read(file, "header/header")
        header <- header[grepl("#samheader: @SQ", header)]
        header <- sub("#samheader: ", "", header)
    }
    ## loading data
    getPath <- function(root, ...){
        paste(root, ..., sep="/")
    }
    createReadsName <- function(ids, width=6, prefix="r"){
        paste0(prefix, formatC(ids, width=width, flag="0"))
    }
    filterByPeak <- function(pos, strand, chr, peaks){
        gi <- GInteractions(anchor1=GRanges(chr, IRanges(pos[, 1], width=150)),
                            anchor2=GRanges(chr, IRanges(pos[, 2], width=150)))
        ol <- findOverlaps(gi, peaks)
        keep <- sort(unique(queryHits(ol)))
        list(pos=pos[keep, , drop=FALSE], strand=strand[keep, , drop=FALSE])
    }
    createAlginment <- function(file, p, idx, width, peaks){
        chr_ <- strsplit(p, "/")[[1]]
        chr_id <- which(chr_=="data")[1]
        chr1 <- chr_[chr_id+1]
        chr2 <- chr_[chr_id+2]
        if(chr1!=chr2){
            return(NULL)
        }
        pos <- h5read(file, getPath(p, "position"))
        strand <- h5read(file, getPath(p, "strand"))
        fil <- filterByPeak(pos, strand, chr1, peaks)
        pos <- fil[["pos"]]
        strand <- fil[["strand"]]
        if(nrow(pos)){
            name <- createReadsName(idx+seq.int(nrow(pos)), width=width)
            flag <- ifelse(strand[, 1]=="-", 16, 0)
            posL <- rowMins(pos)
            isize <- abs(pos[, 1] - pos[, 2])
            cigar <- paste0("100M", isize, "N100M")
            aln <- paste(name, flag, chr1, posL, "50", cigar, "*", 0, isize+201, "*", "*", sep = "\\t")
        }else{
            return(NULL)
        }
    }
    exportBamFile <- function(file, peaks){
        con <- sub(".h5\$", "", file)
        sam_path <- sub("h5\$", "sam", file, ignore.case = TRUE)
        if(sam_path==con){
            sam_path <- paste0(con, ".sam")
        }
        sam_con <- file(sam_path, "w")
        on.exit(close(sam_con))
        ## write header
        header <- getHeader(file)
        writeLines(header, sam_con)
        ## write data
        total <- h5read(file, "header/total")
        total_n <- nchar(total)
        h5content <- h5ls(file)
        h5content <- h5content[, "group"]
        h5content <- h5content[grepl("data.*\\\\d+_\\\\d+", h5content)]
        idx <- 0
        for(p in h5content){
            aln <- createAlginment(file, p, idx, total_n, peaks)
            if(length(aln)){
                writeLines(aln, sam_con)
                idx <- idx + length(aln)
            }
        }
        close(sam_con)
        h5closeAll()
        on.exit()
        si <- do.call(rbind, strsplit(header, "\\\\t"))
        si <- as.numeric(sub("LN:", "", si[, 3]))
        si <- si[!is.na(si)]
        if(length(si)){
            si <- any(si>536870912)
        }else{
            si <- TRUE
        }
        if(si){
            bam <- asBam(sam_path, con,
                        overwrite = TRUE, indexDestination = FALSE)
        }else{
            bam <- asBam(sam_path, con,
                        overwrite = TRUE, indexDestination = TRUE)
        }
        unlink(sam_path)
        invisible(bam)
    }

    peaks <- read.csv(peaks)
    peaks <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)),
                                        GRanges(chr2, IRanges(start2, end2))))
    # output
    null <- lapply(pairs, exportBamFile, peaks=peaks)
    """
}
