// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIR2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(peak), path(pairs)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai")    , emit: bam
    path "versions.yml"                                  , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Oct. 2021 convert pairs.gz to bam file for visualization
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("Rsamtools", "InteractionSet")
    versions <- c("${getProcessName(task.process)}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    peaks <- "$peak"
    pairs <- dir(".", "unselected.pairs.gz")
    cons<- sub("unselected.pairs.gz", "bam", pairs)

    ## scan file header
    scanPairHeader <- function(file){
        con <- gzcon(file(file, open="rb"))
        on.exit(close(con))
        header <- NULL
        while(length(line <- readLines(con, n = 1))){
            if(grepl("^#", line)){
                header <- c(header, line)
            }else{
                return(header)
            }
        }
    }

    ## keep the reads only in interactions
    dataInPeak <- function(data, peak){
        dgi <- with(data, GInteractions(GRanges(V2, IRanges(V3, width=150)),
                                        GRanges(V4, IRanges(V5, width=150))))
        data[countOverlaps(dgi, peak, use.region="same")>0, , drop=FALSE]
    }

    exportBamFile <- function(header, data, con){
        sam_path <- sub("bam\$", "sam", con, ignore.case = TRUE)
        if(sam_path==con){
            sam_path <- paste0(con, ".sam")
        }
        sam_con <- file(sam_path, "w")
        on.exit(close(sam_con))
        data <- data[data[, 2]==data[, 4], , drop=FALSE]
        colnames(data) <- c("name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")
        data\$flag <- ifelse(data\$strand1=="-", 16, 0)
        data\$pos <- rowMins(as.matrix(data[, c("pos1", "pos2")]))
        data\$isize <- abs(data\$pos1 - data\$pos2)
        data\$cigar <- paste0("100M", data\$isize, "N100M")
        aln <- paste(data\$name, data\$flag, data\$chr1, data\$pos, "50", data\$cigar, "*", 0, data\$isize+201, "*", "*", sep = "\t")

        writeLines(header, sam_con)
        writeLines(aln, sam_con)
        close(sam_con)
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
            bam <- asBam(sam_path, sub(".bam\$", "", con, ignore.case = TRUE),
                        overwrite = TRUE, indexDestination = FALSE)
        }else{
            bam <- asBam(sam_path, sub(".bam\$", "", con, ignore.case = TRUE),
                        overwrite = TRUE, indexDestination = TRUE)
        }
        unlink(sam_path)
        invisible(bam)
    }

    peaks <- read.csv(peaks)
    peaks <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)),
                                        GRanges(chr2, IRanges(start2, end2))))
    # output
    null <- mapply(function(file, con){
        header <- scanPairHeader(file)
        header <- header[grepl("#samheader: @SQ", header)]
        header <- sub("#samheader: ", "", header)

        data <- read.table(file,
                        colClasses=c("character", "character",
                                    "integer", "character",
                                    "integer", "character",
                                    "character",
                                    rep("NULL", 7)))
        data <- data[!duplicated(data[, -1]), , drop=FALSE]
        data <- dataInPeak(data, peaks)
        exportBamFile(header, data, con)
    }, pairs, cons)
    """
}
