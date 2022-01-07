process PREPARE_COUNTS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(r2peak, stageAs:"R2peak/*"), path(r1peak, stageAs: "R1peak/*"), path(distalpair, stageAs: "pairs/*"), path(mappability), path(fasta), path(cut)

    output:
    tuple val(meta), path("*.csv"), emit: counts
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: "1e9"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 count reads for peak filtering
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("rtracklayer", "InteractionSet", "Biostrings", "Rsamtools", "rhdf5", "BiocParallel")
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
    OUTPUT <- "counts.${meta.id}.csv"
    NCORE <- as.numeric("$task.cpus")
    peak_pair_block <- ifelse(!is.na(as.numeric($args)), as.numeric($args)[1], 1e9)
    FASTA <- "$fasta"
    CUT <- "$cut"
    MAPPABILITY <- "$mappability"
    pattern <- "h5" ## h5 is postfix of output of pairtools pairs2hdf5
    pairs <- dir("pairs", paste0(pattern, "\$"), full.names=TRUE)
    names(pairs) <- sub(paste0("\\\\.", pattern), "", basename(pairs))
    R1PEAK <- import("$r1peak")
    R2PEAK <- import("$r2peak")
    mcols(R1PEAK) <- NULL
    mcols(R2PEAK) <- NULL
    ## split by chromsome
    R1PEAK <- split(R1PEAK, seqnames(R1PEAK))
    R2PEAK <- split(R2PEAK, seqnames(R2PEAK))
    R1PEAK <- R1PEAK[lengths(R1PEAK)>0]
    R2PEAK <- R2PEAK[lengths(R2PEAK)>0]
    chromosomes <- intersect(names(R1PEAK), names(R2PEAK))
    chromosomes <- chromosomes[!grepl("_", chromosomes)]
    chromosomes <- chromosomes[!grepl("M", chromosomes)] ## remove chrM/chrMT
    if(length(chromosomes)==0){
        stop("no valid data in same chromosome.")
    }

    ## loading data
    readPairs <- function(pair, chrom1, chrom2){
        h5content <- rhdf5::h5ls(pair)
        h5content <- h5content[, "group"]
        h5content <- h5content[grepl("data.*\\\\d+_\\\\d+", h5content)]
        h5content <- unique(h5content)
        n <- h5content[grepl(paste0("data.", chrom1, ".", chrom2), h5content)]
        n <- paste(n, "position", sep="/")
        inf <- rhdf5::H5Fopen(pair, flags="H5F_ACC_RDONLY")
        on.exit({rhdf5::H5Fclose(inf)})
        pc <- lapply(n, function(.ele){
            if(rhdf5::H5Lexists(inf, .ele)){
                rhdf5::h5read(inf, .ele)
            }
        })
        rhdf5::H5Fclose(inf)
        rhdf5::h5closeAll()
        on.exit()
        pc <- do.call(rbind, pc)
    }

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
    gis <- NULL

    gc(reset=TRUE)
    param <- SnowParam(workers = NCORE, progressbar = TRUE, type = "SOCK")
    for(chrom1 in chromosomes){
        for(chrom2 in chromosomes){
            message("working on ", chrom1, " and ", chrom2)
            r1peak <- R1PEAK[[chrom1]]
            r2peak <- R2PEAK[[chrom2]]
            message("read reads")
            reads <- bplapply(pairs, readPairs, chrom1=chrom1, chrom2=chrom2, BPPARAM = param)
            h5closeAll()
            reads <- do.call(rbind, c(reads, make.row.names = FALSE))
            if(length(reads)){
                reads <- GInteractions(GRanges(chrom1, IRanges(reads[, 1], width=150)),
                                        GRanges(chrom2, IRanges(reads[, 2], width=150)))
                peak_pairs <- expand.grid(seq_along(r1peak), seq_along(r2peak))
                peak_pairs <- split(peak_pairs,
                                    rep(seq.int(ceiling(nrow(peak_pairs)/peak_pair_block)),
                                        each=peak_pair_block)[seq.int(nrow(peak_pairs))])
                message("count reads")
                gi <- bplapply(peak_pairs, FUN=function(peak_pair, reads, r1peak, r2peak){
                    .gi <- InteractionSet::GInteractions(r1peak[peak_pair[, 1]], r2peak[peak_pair[, 2]])
                    S4Vectors::mcols(.gi)[, "count"] <- InteractionSet::countOverlaps(.gi, reads, use.region="both")
                    S4Vectors::mcols(.gi)[, "shortCount"] <- GenomicRanges::countOverlaps(S4Vectors::second(.gi), S4Vectors::second(reads))
                    .gi[S4Vectors::mcols(.gi)[, "count"]>0 & S4Vectors::mcols(.gi)[, "shortCount"]>0]
                }, reads=reads, r1peak=r1peak, r2peak=r2peak, BPPARAM = param)
                gi <- Reduce(c, gi)
                if(length(gi)>0){
                    gis <- c(gis, gi)
                }
                rm(peak_pairs)
                gc(reset=TRUE)
            }
        }
    }
    gis <- do.call(c, gis)
    ### load gc content
    fa <- FaFile(file=FASTA)
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
