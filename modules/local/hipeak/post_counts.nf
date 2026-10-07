process POST_COUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(counts), path(mappability), path(fasta), path(cut)

    output:
    tuple val(meta), path("*.csv"), emit: counts
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 count reads for peak filtering
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("rtracklayer", "InteractionSet", "Biostrings", "Rsamtools")
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
    ## make_option(c("-a", "--r1peak"), type="character", default=NULL, help="filename of r1 peak", metavar="string")
    ## make_option(c("-b", "--r2peak"), type="character", default=NULL, help="filename of r2 peak", metavar="string")
    ## make_option(c("-x", "--restrict"), type="character", default=NULL, help="filename of restrict cut", metavar="string")
    ## make_option(c("-p", "--pairs"), type="character", default=NULL, help="folder of valid distal pairs", metavar="string")
    ## make_option(c("-m", "--mappability"), type="character", default=NULL, help="mappability file", metavar="string")
    ## make_option(c("-o", "--output"), type="character", default="counts.csv", help="output folder", metavar="string")
    ## make_option(c("-f", "--fasta"), type="character", default=NULL, help="genome fasta file", metavar="string")
    ## make_option(c("-1", "--chrom1"), type="character", default=NULL, help="chromosome1", metavar="string")
    ## make_option(c("-2", "--chrom2"), type="character", default=NULL, help="chromosome1", metavar="string")
    OUTPUT <- "counts.${meta.id}.csv"
    FASTA <- "$fasta"
    CUT <- "$cut"
    MAPPABILITY <- "$mappability"

    getMscore <- function(gr){
        gr1 <- split(unique(gr), as.character(seqnames(unique(gr))))
        gr0 <- lapply(gr1, function(.ele) .ele[seq.int(min(50, length(.ele)))])
        available_chr <- vapply(gr0, FUN=function(chr){
            out <- try(ms <- import(MAPPABILITY, which=chr))
            if(inherits(out, "try-error")) return(FALSE)
            return(length(ms)>0)
        }, FUN.VALUE=logical(1))
        available_chr <- names(available_chr)[available_chr]
        mscore <- import(MAPPABILITY, which=gr[seqnames(gr) %in% available_chr], as="RleList")
        chr <- intersect(names(mscore), names(gr1))
        vw <- Views(mscore[chr], gr1[chr])
        sc <- viewMeans(vw)
        vs <- ranges(vw)
        vs <- as(vs, "GRanges")
        vs\$score <- unlist(sc, use.names=FALSE)
        ol <- findOverlaps(gr, vs, type="equal")
        ol <- ol[!duplicated(queryHits(ol))]
        score <- vs[subjectHits(ol)]\$score[match(seq_along(gr), queryHits(ol))]
        score[is.na(score)] <- 0
        score
    }

    ### load counts
    gis_rds <- dir(".", pattern="counts.${meta.id}.*.rds")
    gis <- lapply(gis_rds, readRDS)
    gis <- gis[lengths(gis)>0]
    gis <- do.call(c, gis)
    stopifnot(is(gis, "GInteractions"))
    ### load gc content
    fa_idx <- indexFa(file=FASTA)
    fa <- FaFile(file=FASTA, index=fa_idx)
    seqinfo(gis) <- seqinfo(fa)[seqlevels(gis)]
    gis <- trim(gis)
    gc1 <- letterFrequency(getSeq(fa, first(gis)), letters="CG", as.prob=TRUE)
    gc2 <- letterFrequency(getSeq(fa, second(gis)), letters="CG", as.prob=TRUE)
    gis\$gc <- gc1 * gc2 + 1e-9
    ### load enzyme cut number in R1
    cut <- import(CUT)
    start(cut) <- end(cut)
    gis\$cut <- countOverlaps(first(gis), cut)+0.1
    ### load mapping score
    m1 <- getMscore(first(gis))
    m2 <- getMscore(second(gis))
    gis\$mappability <- m1*m2 + 1e-6
    ### get distance of the anchors
    gis\$dist <- distance(first(gis), second(gis))+1
    gis\$dist[is.na(gis\$dist)] <- max(gis\$dist, na.rm=TRUE)*100 ##trans interactions
    gis\$length <- width(first(gis))*width(second(gis))

    mm <- log(as.data.frame(mcols(gis)[, c("length", "cut", "gc", "mappability", "dist", "shortCount")]))
    mm <- cbind(as.data.frame(first(gis)), as.data.frame(second(gis)), gis\$count, mm)
    colnames(mm) <- c("chr1", "start1", "end1", "width1", "strand1",
                    "chr2", "start2", "end2", 'width2', "strand2",
                    "count", "logl", "logn", "loggc", "logm", "logdist", 'logShortCount')
    mm\$strand1 <- NULL
    mm\$strand2 <- NULL
    write.csv(mm, OUTPUT, row.names=FALSE)
    """
}
