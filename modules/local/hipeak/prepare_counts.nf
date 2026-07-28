process PREPARE_COUNTS {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(r2peak, stageAs:"R2peak/*"), path(r1peak, stageAs: "R1peak/*"), path(distalpair, stageAs: "pairs/*"), val(chrom1)

    output:
    tuple val(meta), path("*.rds"), optional: true, emit: counts
    path "versions.yml"                           , emit: versions

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
    pkgs <- c("rtracklayer", "InteractionSet", "rhdf5", "BiocParallel")
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
    parse_args <- function(options, args){
        out <- lapply(options, function(.ele){
            if(any(.ele[-3] %in% args)){
                if(.ele[3]=="logical"){
                    TRUE
                }else{
                    id <- which(args %in% .ele[-3])[1]
                    x <- args[id+1]
                    mode(x) <- .ele[3]
                    x
                }
            }
        })
    }
    option_list <- list("peak_pair_block"=c("--peak_pair_block", "-b", "integer"),
                        "snow_type"=c("--snow_type", "-t", "character"))
    opt <- parse_args(option_list, strsplit("$args", "\\\\s+")[[1]])
    CHROM1 <- "$chrom1"
    OUTPUT <- "counts.${meta.id}.${chrom1}.rds"
    NCORE <- as.numeric("$task.cpus")
    SNOW_TYPE <- "SOCK"
    COUNTS_FILTER <- 1
    peak_pair_block <- 1e9
    if(!is.null(opt\$peak_pair_block)){
        peak_pair_block <- opt\$peak_pair_block
    }
    if(!is.null(opt\$snow_type)){
        SNOW_TYPE <- opt\$snow_type
    }
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
    chromosomes <- chromosomes[!grepl("EBV", chromosomes)]
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

    ### load counts
    gis <- NULL

    if(CHROM1 %in% chromosomes){
        gc(reset=TRUE)
        if(SNOW_TYPE=="FORK"){
            param <- MulticoreParam(workers = NCORE, progressbar = TRUE)
        }else{
            param <- SnowParam(workers = NCORE, progressbar = TRUE, type = SNOW_TYPE)
        }
        chrom1 <- CHROM1
        parallel <- TRUE
        for(chrom2 in chromosomes){
            message("working on ", chrom1, " and ", chrom2, " from ", Sys.time())
            r1peak <- R1PEAK[[chrom1]]
            r2peak <- R2PEAK[[chrom2]]
            message("read reads")
            if(parallel){
                try_res <- try({reads <- bplapply(pairs, readPairs, chrom1=chrom1, chrom2=chrom2, BPPARAM = param)})
                if(inherits(try_res, "try-error")){
                    parallel <- FALSE
                }
            }
            if(!parallel){
                reads <- lapply(pairs, readPairs, chrom1=chrom1, chrom2=chrom2)
            }
            h5closeAll()
            reads <- do.call(rbind, c(reads, make.row.names = FALSE))
            if(length(reads) && length(r1peak) && length(r2peak)){
                reads <- GInteractions(GRanges(chrom1, IRanges(reads[, 1], width=150)),
                                        GRanges(chrom2, IRanges(reads[, 2], width=150)))
                ## count twice,
                ## first, merge the r1peak with gap 5k and filter the peak with 0 counts
                ## second, count for filtered r1peak
                r1peak_s <- reduce(r1peak, min.gapwidth=5000, with.revmap=TRUE)
                countFUN <- function(peak_pair, reads, r1peak, r2peak){
                    r1peak\$revmap <- NULL
                    .gi <- InteractionSet::GInteractions(r1peak[peak_pair[, 1]], r2peak[peak_pair[, 2]],
                                                        p1=peak_pair[, 1], p2=peak_pair[, 2])
                    reads <- IRanges::subsetByOverlaps(reads, InteractionSet::regions(.gi))
                    ## remove the interactions with distance smaller than 1K
                    .dist <- IRanges::distance(first(.gi), second(.gi))
                    .dist[is.na(.dist)] <- 3e9
                    S4Vectors::mcols(.gi)[, "count"] <- InteractionSet::countOverlaps(.gi, reads, use.region="both")
                    S4Vectors::mcols(.gi)[, "shortCount"] <- GenomicRanges::countOverlaps(S4Vectors::second(.gi), S4Vectors::second(reads))
                    .gi[S4Vectors::mcols(.gi)[, "count"]>=COUNTS_FILTER & S4Vectors::mcols(.gi)[, "shortCount"]>=COUNTS_FILTER & .dist>1000]
                }
                countFUNbyPairs <- function(r1peak, r2peak, peak_pairs, reads, parallel){
                    peak_pairs_group <- ceiling(nrow(peak_pairs)/peak_pair_block)
                    if(peak_pairs_group>1){
                        peak_pairs_group <- rep(seq.int(peak_pairs_group), each= peak_pair_block)
                        peak_pairs <- split(peak_pairs,
                                            peak_pairs_group[seq.int(nrow(peak_pairs))])
                        if(parallel){
                            try_res <- try({gi <- bplapply(peak_pairs, FUN=countFUN, reads=reads, r1peak=r1peak, r2peak=r2peak, BPPARAM = param)})
                            if(inherits(try_res, "try-error")){
                                parallel <- FALSE
                            }
                        }
                        if(!parallel){
                            gi <- lapply(peak_pairs, FUN=countFUN, reads=reads, r1peak=r1peak, r2peak=r2peak)
                        }
                        gi <- Reduce(c, gi)
                    }else{
                        gi <- countFUN(peak_pairs, reads, r1peak, r2peak)
                    }
                    gi
                }
                message("count reads")
                peak_pairs <- expand.grid(seq_along(r1peak_s), seq_along(r2peak))
                gi <- countFUNbyPairs(r1peak_s, r2peak, peak_pairs, reads, parallel)
                peak_pairs <- mcols(r1peak_s)[mcols(gi)[, "p1"], "revmap"]
                peak_pairs <- data.frame(p1=unlist(peak_pairs),
                                        p2=rep(mcols(gi)[, "p2"], lengths(peak_pairs)))
                rm(gi)
                gi <- countFUNbyPairs(r1peak, r2peak, peak_pairs, reads, parallel)
                if(length(gi)>0){
                    gis <- c(gis, gi)
                }
                rm(peak_pairs, gi)
                gc(reset=TRUE)
            }
        }
    }
    if(is.list(gis)) gis <- do.call(c, gis)
    saveRDS(gis, OUTPUT)
    """
}
