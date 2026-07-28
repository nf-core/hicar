process BIOC_ENZYMECUT {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(bin_size), val(site), path(fasta)
    val enzyme

    output:
    tuple val(bin_size), path('*.cut')        , emit: cut
    path "versions.yml"                       , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 for trackViewer parser
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################

    pkgs <- c("Biostrings", "GenomicRanges", "BSgenome")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    ## enzyme-cut
    fa <- readDNAStringSet("${fasta}")
    names(fa) <- sub("\\\\s+.*\$", "", names(fa))
    site <- strsplit("$site", "\\\\s+")[[1]]
    pattern <- DNAString(site[1])
    pos <- as.numeric(site[2])
    hits <- vmatchPattern(pattern, fa)

    #hits <- as(hits, "GRanges")
    #seqinfo(hits) <- seqinfo(fa)[seqlevels(hits)]
    gaps <- mapply(hits, seqlengths(fa)[names(hits)], FUN=function(.hits, .seql){
        hits_l <- shift(.hits, shift=pos)
        width(hits_l) <- 1
        hits_r <- .hits
        start(hits_r) <- end(hits_r) - pos
        width(hits_r) <- 1

        frag_l_pos <- start(hits_l)
        frag_r_pos <- start(hits_r)
        frag_l_ir <- IRanges(frag_l_pos, c(frag_l_pos[-1], .seql))
        mcols(frag_l_ir)\$strand <- "+"
        frag_r_ir <- IRanges(c(1, frag_r_pos[-length(frag_r_pos)]), frag_r_pos)
        mcols(frag_r_ir)\$strand <- "-"
        ir <- sort(c(frag_r_ir, frag_l_ir))
        mcols(ir)\$idx <- seq_along(ir)
        ir
    }, SIMPLIFY = FALSE)

    gaps <- mapply(gaps, names(gaps), FUN=function(ir, seq){
        strand <- mcols(ir)\$strand
        mcols(ir)\$strand <- NULL
        GRanges(rep(seq, length(ir)), ir, strand=strand)
    })
    gaps <- unlist(GRangesList(gaps))
    seqinfo(gaps) <- seqinfo(fa)[seqlevels(gaps)]
    gaps_promters <- promoters(gaps, upstream=0, downstream=200)
    gaps_promters <- trim(gaps_promters)
    gaps_promters <- getSeq(fa, gaps_promters)
    gaps\$gc <- letterFrequency(gaps_promters, letters="GC", OR="|", as.prob=TRUE)

    out <- as.data.frame(unname(gaps))
    out\$pos <- ifelse(out\$strand=="-", end(gaps), start(gaps))
    out <- out[, c("idx", "strand", "seqnames", "pos", "width", "G.C")]
    write.table(out, paste0("$enzyme", ".cut"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    """
}
