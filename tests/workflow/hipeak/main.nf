#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// run test: PROFILE=docker pytest --tag callhipeak --symlink --kwdof
// run simulate: PROFILE=docker nextflow run path/to/nf-core-hicar/tests/workflow/hipeak/ -entry test_call_hi_peak -c path/to/nf-core-hicar/nextflow.config,path/to/nf-core-hicar/tests/config/nextflow.config --reads 1e6

include { GUNZIP
    } from '../../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP2
    } from '../../../modules/nf-core/modules/gunzip/main'
include { CHROMSIZES
    } from '../../../modules/local/genome/chromsizes'
include { COOLER_DIGEST
    } from '../../../modules/nf-core/modules/cooler/digest/main'
include { BIOC_PAIRS2HDF5
    } from '../../../modules/local/bioc/pairs2hdf5'
include { ATAC_PEAK
    } from '../../../subworkflows/local/callatacpeak'
include { R1_PEAK
    } from '../../../subworkflows/local/calldistalpeak'
include { HI_PEAK
    } from '../../../subworkflows/local/hipeak'

process CREATE_PAIRS {
    publishDir "${params.outdir}/pairs",
        mode: "copy"
    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    path chromsizes
    path gtf
    path digest_genome_bed
    val reads
    val seed

    output:
    path "$peak1"       , emit: peak1
    path "$peak2"       , emit: peak2
    path "$pairs"       , emit: pairs
    path "$distalpairs" , emit: distalpairs
    path "$interactions", emit: interactions

    script:
    peak1="peak1.bed"
    peak2="peak2.bed"
    pairs="test.unselected.pairs.gz"
    distalpairs="distal.pairs.gz"
    cut="CviQI.bed"
    interactions="interaction.bedpe"
    """
    #!/usr/bin/env Rscript
    library(GenomicRanges)
    library(rtracklayer)
    library(GenomicFeatures)
    library(InteractionSet)
    options(scipen=10)
    set.seed($seed)
    peak_num_1 <- 10000      # fake R1 events
    peak_num_2 <- 1000       # fake R2 events
    real_r1 <- 1000          # real R1 events
    real_r2 <- 100           # real R2 events
    real_conn_num  <- 1e4    # real connection events
    peak_reads_num <- $reads # total signal reads
    background_lambda <- .5  # background lambda
    real_connection <- expand.grid(sample.int(peak_num_1, real_r1), sample.int(peak_num_2, real_r1))
    real_connection <- real_connection[sample.int(nrow(real_connection), real_conn_num), ]
    rtnorm <- function(n, mean = 0, sd = 1, min = 0, max = 1) { ## help function
        bounds <- pnorm(c(min, max), mean, sd)
        u <- runif(n, bounds[1], bounds[2])
        round(qnorm(u, mean, sd))
    }
    ## R1 peak candidates
    r1s <- import("$digest_genome_bed")
    r1s <- r1s[-1]
    width(r1s) <- 1
    export(r1s, "$cut")
    chromsize <- read.delim("$chromsizes", header=FALSE)
    chromosomes <- chromsize[, 1]
    seqlengths <- chromsize[, 2]
    names(seqlengths) <- chromosomes
    tileGenome <- tileGenome(seqlengths, tilewidth=500)
    tileGenome <- unlist(tileGenome)
    tileGenome\$count <- countOverlaps(tileGenome, r1s)
    tileGenome <- tileGenome[tileGenome\$count>0]
    stopifnot(length(tileGenome)>peak_num_1*2)
    peaks1 <- sample(seq_along(tileGenome), peak_num_1*2, replace=FALSE)
    peaks1 <- tileGenome[peaks1]
    peaks1 <- reduce(peaks1)
    peaks1 <- peaks1[sample.int(length(peaks1), peak_num_1)]

    ## R2 peak candidates from promoter
    txdb <- makeTxDbFromGFF("$gtf")
    gene <- genes(txdb)
    pro <- promoters(gene, upstream=4000, downstream=1000)
    pro <- reduce(pro)
    stopifnot(length(pro)>peak_num_2)
    peaks2 <- sample.int(length(pro), peak_num_2, replace=FALSE)
    peaks2 <- pro[peaks2]
    peaks2 <- shift(peaks2, shift=2000)
    w <- rtnorm(peak_num_2, mean=800, sd=300, min=350, max=5000)
    width(peaks2) <- w

    ## real connections
    real_gi <-GInteractions(anchor1=peaks1[real_connection[, 1]],
                            anchor2=peaks2[real_connection[, 2]])
    ## generate R1 background reads
    ol <- findOverlaps(peaks1, r1s)
    q1 <- unique(subjectHits(ol))
    bg_ids <- rpois(length(q1), lambda=background_lambda) #background for r1
    bg_ids <- rep(q1, bg_ids)
    r1_background <- r1s[bg_ids]
    ## generate R2 background reads
    q2 <- tile(peaks2, width=1)
    q2_bkg <- unlist(q2)
    q2_bkg <- q2_bkg[sample.int(length(q2_bkg), length(r1s))]
    r2_background <- q2_bkg[bg_ids]
    ## generate ATAC reads, r1/r2 distance < 1e4
    ol <- findOverlaps(peaks2, r1s, maxgap=1e3)
    q1 <- split(subjectHits(ol), queryHits(ol))
    keep <- unique(queryHits(ol))
    mn <- peak_reads_num*.6/length(peaks2)
    atac_ids <- rtnorm(length(keep), mean = mn, sd = mn/10, min=mn/2, max=mn*10) # 60%
    atac_r2_ids <- mapply(lengths(q2[keep]), atac_ids, FUN=sample.int, replace=TRUE, SIMPLIFY=FALSE)
    atac_r2_ids <- c(atac_r2_ids[[1]], mapply(cumsum(lengths(q2[keep][-length(q2[keep])])), atac_r2_ids[-1], FUN=`+`, SIMPLIFY=FALSE))
    r2_atac <- unlist(q2[keep])[unlist(atac_r2_ids)]
    atac_r1_ids <- mapply(lengths(q1), atac_ids, FUN=sample.int, replace=TRUE, SIMPLIFY=FALSE)
    atac_r1_ids <- c(atac_r1_ids[[1]], mapply(cumsum(lengths(q1[-length(q1)])), atac_r1_ids[-1], FUN=`+`, SIMPLIFY=FALSE))
    r1_atac <- r1s[unlist(q1)[unlist(atac_r1_ids)]]

    ## generate real reads
    mn <- peak_reads_num*.4/real_conn_num
    mcols(real_gi)[, "counts"] <- rtnorm(real_conn_num, mean=mn, sd=mn/2, min=mn/10, max=mn*10)
    ol <- findOverlaps(first(real_gi), r1s)
    ols <- split(subjectHits(ol), queryHits(ol))
    r1_ids <- mapply(lengths(ols), mcols(real_gi)[, "counts"], FUN=sample, replace=TRUE, SIMPLIFY=FALSE)
    r1_ids <- mapply(ols, r1_ids, FUN=function(a, i) a[i], SIMPLIFY=FALSE)
    r1_ids <- unlist(r1_ids)
    r1_signal <- r1s[r1_ids]
    q2 <- tile(second(real_gi), width=50)
    r2_ids <- mapply(lengths(q2), mcols(real_gi)[, "counts"], FUN=sample, replace=TRUE, SIMPLIFY=FALSE)
    r2_ids <- c(r2_ids[[1]], mapply(cumsum(lengths(q2[-length(q2)])), r2_ids[-1], FUN=`+`, SIMPLIFY=FALSE))
    r2_signal <- unlist(q2)[unlist(r2_ids)]
    r2_signal <- shift(r2_signal, shift=sample.int(50, length(r2_signal), replace=TRUE))
    width(r2_signal) <- 1
    reads1 <- c(r1_background, r1_atac, r1_signal)
    reads2 <- c(r2_background, r2_atac, r2_signal)
    ## sort the reads by r1
    ord <- order(reads1, reads2)
    reads1 <- reads1[ord]
    reads2 <- reads2[ord]

    N <- length(reads1)
    readID <- paste0("r", formatC(seq.int(N), width=nchar(as.character(N)), flag="0"))
    strand1 <- sample(c("+", "-"), N, replace=TRUE)
    strand2 <- ifelse(strand1=="-", "+", "-")

    # write peaks1.bed
    export(peaks1, "$peak1")
    # write peaks2.bed
    export(peaks2, "$peak2")
    # write pairs
    header <- c("## pairs format v1.0.0",
                "#shape: whole matrix",
                "#genome_assembly: unknown",
                paste("#chromosomes:", paste(chromsize[, 1], collapse=" ")),
                paste("#chromsize:", chromsize[, 1], chromsize[, 2]),
                paste0("#samheader: @SQ	SN:", chromsize[, 1],"	LN:", chromsize[, 2]),
                "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type")
    content <- paste(readID,
                    as.character(seqnames(reads1)),
                    start(reads1),
                    as.character(seqnames(reads2)),
                    start(reads2),
                    strand1, strand2, "UU",
                    sep="\\t")
    gz <- gzfile("$pairs", "w")
    writeLines(c(header, content), gz)
    close(gz)
    # write distalpair
    dist <- distance(reads1, reads2)
    content_dist <- content[dist>2000]
    gz2 <- gzfile("$distalpairs", "w")
    writeLines(c(header, content_dist), gz2)
    close(gz2)
    # write real connections
    export(real_gi, "$interactions")
    """
}

process CHECK_PEAKS {
    publishDir "${params.outdir}",
        mode: "copy"
    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }
    input:
    tuple val(meta), path(peak)
    path peak2
    path peak1
    path pairs
    path interactions
    path hdf5

    output:
    path "results.txt", emit: qc
    path "results.csv", emit: res

    script:
    """
    #!/usr/bin/env Rscript
    library(rtracklayer)
    library(InteractionSet)
    library(rhdf5)
    hdf5 <- "$hdf5"
    peaks1 <- import("$peak1")
    peaks2 <- import("$peak2")
    reads <- read.table("$pairs")
    h5content <- h5ls(hdf5)
    h5content <- h5content[, "group"]
    h5content <- h5content[grepl("data.*\\\\d+_\\\\d+", h5content)]
    h5content <- unique(h5content)
    h5reads_pos <- lapply(paste0(h5content, "/position"), h5read, file=hdf5)
    h5reads_pos <- do.call(rbind, h5reads_pos)
    stopifnot("There are bugs in creating pairs2hdf5.nf"=identical(nrow(h5reads_pos), nrow(reads)))
    peaks <- read.table("$peak", header=TRUE)
    gi <- with(reads, GInteractions(GRanges(V2, IRanges(V3, width=150)), GRanges(V4, IRanges(V5, width=150))))
    peaks_gi <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)), GRanges(chr2, IRanges(start2, end2))))
    peaks1\$count <- countOverlaps(peaks1, gi)
    peaks2\$count <- countOverlaps(peaks2, gi)
    peaks_gi\$count <- countOverlaps(peaks_gi, gi, use.region="both")
    pct <- sum(peaks_gi\$count==peaks\$count)/length(peaks_gi)
    stopifnot("There are bugs in counting"=pct>.99)
    real_gi <- import("$interactions", format="bedpe")
    detected_r <- subsetByOverlaps(real_gi, peaks_gi)
    detected_p <- subsetByOverlaps(peaks_gi, real_gi)
    P <- length(real_gi)                             ## condition positive
    N <- length(peaks1)*length(peaks2) - P           ## condition negative
    TP <- length(detected_r)                         ## True positive
    FP <- length(peaks_gi) - length(detected_p)      ## False positive
    sensitivity <- TP/P                              ## sensitivity
    FDR <- FP/length(peaks_gi)                       ## false discovery rate
    write.csv(c("Total reads"=nrow(reads),
                "Total true connections"=P,
                "True positive"=TP,
                "False positive"=FP,
                "sensitivity"=sensitivity),
                "results.csv")
    writeLines(ifelse(FDR<0.1, "YES", "NO"), "results.txt")
    """
}

workflow test_call_hi_peak {
    fasta           = GUNZIP([[id:"test"],file(params.test_data.fasta, checkIfExists: true)]).gunzip.map{it[1]}
    chromsizes      = CHROMSIZES ( fasta ).sizes
    macs_gsize      = params.test_data.macs_gsize
    pval            = params.r1_pval_thresh
    gtf             = file(params.test_data.gtf, checkIfExists: true)
    mappability     = file(params.test_data.mappability, checkIfExists: true)
    digest_genome_bed = COOLER_DIGEST (
        fasta,
        chromsizes,
        params.enzyme
    ).bed
    CREATE_PAIRS(chromsizes, gtf, digest_genome_bed, params.reads, params.seed)
    ATAC_PEAK(
        GUNZIP2(CREATE_PAIRS.out.pairs.map{[[id:"test", group:"gp1"], it]}).gunzip,
        chromsizes,
        macs_gsize,
        gtf,
        "HiCAR",
        null,
        Channel.empty()
    )
    R1_PEAK(
        CREATE_PAIRS.out.distalpairs.map{[[id:"test", group:"gp1"], it]},
        chromsizes,
        digest_genome_bed,
        gtf,
        pval
    )
    BIOC_PAIRS2HDF5(
        CREATE_PAIRS.out.pairs.map{[[id:"test", group:"gp1"], it]},
        chromsizes
    )
    hdf5 = BIOC_PAIRS2HDF5.out.hdf5.map{[it[1]]}
    validpair  = ATAC_PEAK.out.mergedpeak
                    .combine( R1_PEAK.out.mergedpeak )
                    .combine( hdf5 )
                    .map{[[id:"test"], it[0], it[1], it[2]]}
    HI_PEAK (
        validpair,
        chromsizes,
        Channel.fromPath(params.test_data.gtf),
        fasta,
        digest_genome_bed,
        Channel.fromPath(params.test_data.mappability),
        true, true
    )
    CHECK_PEAKS (
        HI_PEAK.out.peak,
        CREATE_PAIRS.out.peak2,
        CREATE_PAIRS.out.peak1,
        CREATE_PAIRS.out.pairs,
        CREATE_PAIRS.out.interactions,
        hdf5
    )
}
