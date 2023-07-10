process HICDCPLUS_FEATURES {
    tag "$bin_size"
    label 'process_high'
    label 'error_ignore'

    conda "bioconda::bioconductor-hicdcplus=1.2.1"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-hicdcplus:1.2.1--r41h619a076_0' :
        'biocontainers/bioconductor-hicdcplus:1.2.1--r41h619a076_0' }"

    input:
    tuple val(bin_size), val(site), path(fasta), path(chrom_sizes), path(mappability)

    output:
    tuple val(bin_size), path("*_bintolen.txt.gz")            , emit: features
    path "versions.yml"                                       , emit: versions

    script:
    prefix   = task.ext.prefix ?: "hicdcplus_bin${bin_size}"
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on July 13, 2022 for HiC-DC+ to prepare the features
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################

    pkgs <- c("HiCDCPlus", "dplyr", "BSgenome", "GenomicRanges", "BiocGenerics",
                "Biostrings", "rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    ## arguments
    FASTA <- "$fasta"
    CHROMSIZE_FILE <- "$chrom_sizes"
    BIN <- as.numeric("$bin_size")
    SITES <- "$site"
    MAPPABILITY <- "$mappability"
    bintolenoutput <- "${prefix}_bintolen.txt.gz"
    sig <- strsplit(SITES, "\\\\s+")[[1]][1]
    binsize <- as.numeric(BIN)
    chrom_sizes <- read.delim(CHROMSIZE_FILE, header=FALSE)
    chrs <- chrom_sizes[, 1]
    chrs <- chrs[!grepl("_", chrs)]
    chrs <- chrs[!grepl("EBV|M|Y", chrs)]
    chrom_sizes <- chrom_sizes[match(chrs, chrom_sizes[, 1]), , drop=FALSE]
    chrom_sizes <- chrom_sizes[, 2]
    names(chrom_sizes) <- chrs

    bin_type <- "Bins-uniform"
    feature_type <- "RE-based"
    wg_file <- MAPPABILITY
    ## handle args
    args <- strsplit("${args}", "\\\\s+")[[1]]
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
    option_list <- list("bin_type"=c("--bin_type", "-b", "charcter"),
                        "feature_type"=c("--feature_type", "-f", "character"),
                        "wg_file"=c("--wg_file", "-w", "character"))
    opt <- parse_args(option_list, args)
    if(!is.null(opt\$bin_type)){
        bin_type <- opt\$bin_type
    }
    if(!is.null(opt\$feature_type)){
        feature_type <- opt\$feature_type
    }
    if(!is.null(opt\$wg_file)){
        wg_file <- opt\$wg_file
    }
    genome.chromGR <- GRanges(chrs, IRanges(1, width=chrom_sizes))
    if (!(is.null(wg_file))) {
        wgdata <- import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    }else{
        wgdata <- GRanges()
    }
    wgdata <- subsetByOverlaps(wgdata, genome.chromGR)

    genome <- readDNAStringSet(FASTA)
    names(genome) <- sub("\\\\s+.*\$", "", names(genome))
    stopifnot(all(chrs %in% names(genome)))

    getGC <- function(genome, gr, as.prob=TRUE){
        letterFrequency(getSeq(genome, gr), 'GC', as.prob=as.prob)
    }
    getMapScore <- function(data, regions){
        if(is.null(data)) return(NULL)
        if(length(data)==0) return(NULL)
        if(length(data\$score)==0) return(NULL)
        hits <- GenomicRanges::findOverlaps(data, regions, type = "within")
        DT <- data.frame(queryHits = as.data.frame(hits)[[1]],
                        subjectHits = as.data.frame(hits)[[2]])
        DT <- DT %>% dplyr::mutate(map = data\$score[.data\$queryHits],
                        len = GenomicRanges::end(data)[.data\$queryHits] -
                                GenomicRanges::start(data)[.data\$queryHits] +
                                    1) %>% dplyr::group_by(.data\$subjectHits) %>%
                                        dplyr::summarize(avgmap = stats::weighted.mean(map, w = len))
        map <- rep(0, length(regions))
        map[DT\$subjectHits] <- DT\$avgmap
        return(map)
    }

    if(feature_type == "RE-based"){
        patternobj <- Biostrings::DNAStringSet(sig)
        enzymeCuts <- lapply(genome[chrs], function(x){
            unlist(as(matchPDict(patternobj, x), "CompressedIRangesList"))
        })
        enzymeCuts <- sort(as(as(enzymeCuts, "IRangesList"), "GRanges"))
        seqlengths(enzymeCuts) <- chrom_sizes[seqlevels(enzymeCuts)]
        enzymeCuts <- trim(enzymeCuts)
        RE_sites <- enzymeCuts
        width(RE_sites) <- 1
        RE_sites <- gaps(RE_sites)
        RE_sites <- RE_sites[strand(RE_sites)=="*"]
        ids <- start(RE_sites) != 1
        start(RE_sites[ids]) <- start(RE_sites[ids]) - 1
    }else{
        stopifnot("feature_type=RE-agnostic requires a fixed binsize." = bin_type == "Bins-uniform")
    }

    if (bin_type == "Bins-RE-sites"){
        endL <- endR <- RE_sites
        endL <- trim(resize(endL, 200, fix="start"))
        endR <- trim(resize(endR, 200, fix="end"))
        gcL <- getGC(genome, endL)
        gcR <- getGC(genome, endR)
        gc <- (gcL+gcR)/2
        RE_sites\$gc <- gc
        if(length(wgdata)>0){
            endL <- trim(resize(endL, 500, fix="start"))
            endR <- trim(resize(endR, 500, fix="end"))
            mapL <- getMapScore(regions = endL, data = wgdata)
            mapR <- getMapScore(regions = endR, data = wgdata)
            map <- (mapL + mapR)/2
            if(length(map)==length(RE_sites)){
                RE_sites\$map <- map
            }else{
                RE_sites\$map <- 0
            }
        }else{
            RE_sites\$map <- 0
        }
        bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
        bintolen <- bintolen %>% dplyr::mutate(RE_id = dplyr::row_number())
        if (binsize > 1) {
            bintolen <- bintolen %>%
                dplyr::mutate(
                    binNumb = floor((dplyr::row_number() - 1)/binsize) + 1)
            bintolen <- bintolen %>%
                dplyr::group_by(.data\$seqnames, .data\$binNumb) %>%
                    dplyr::summarize(
                        start = min(start),
                        end = max(end),
                        gc = mean(gc),
                        map = mean(map),
                        RE_id = min(.data\$RE_id)) %>%
                            dplyr::mutate(
                                width = .data\$end - .data\$start + 1) %>%
                                    dplyr::rename(chr = "seqnames") %>%
                                        dplyr::mutate(bins =
                                            paste(.data\$chr,
                                                .data\$start,
                                                .data\$end,
                                                sep = "-"))
        } else {
            bintolen <- bintolen %>% dplyr::rename(chr = "seqnames") %>%
                dplyr::mutate(bins = paste(.data\$chr, .data\$start,
                    .data\$end, sep = "-"))
        }
        bintolen <- bintolen %>% dplyr::select(.data\$bins, .data\$gc,
            .data\$map, .data\$width, .data\$RE_id)
        bintolen <- as.data.frame(bintolen, stringsAsFactors = FALSE)
        rownames(bintolen) <- bintolen\$bins
    }else{
        stopifnot(bin_type=="Bins-uniform")
        binsGR <- tileGenome(chrom_sizes, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
        ids <- start(binsGR) != 1
        start(binsGR[ids]) <- start(binsGR[ids]) - 1
        names(binsGR) <- paste(as.character(seqnames(binsGR)), start(binsGR), end(binsGR), sep = "-")
        binsGR\$map <- binsGR\$len <- binsGR\$gc <- rep(0, length(binsGR))
        if(feature_type == "RE-based"){
            medians <- (start(enzymeCuts) + end(enzymeCuts))/2
            FragmentendsL <- FragmentendsR <- enzymeCuts
            end(FragmentendsL) <- medians - 1
            FragmentendsL <- trim(resize(FragmentendsL, 500, fix="end"))
            start(FragmentendsR) <- medians
            FragmentendsR <- trim(resize(FragmentendsR, 500, fix="start"))
            ends <- c(FragmentendsL, FragmentendsR)
            ol <- findOverlaps(ends, binsGR, type = "within", select = "all")
            LR <- GRanges(names(binsGR[subjectHits(ol)]), ranges(ends[queryHits(ol)]))
            LR <- reduce(LR)
            idx <- as.character(seqnames(LR))
            LR1 <- GRanges(sub("-.*\$", "", idx), ranges(LR), idx=idx)
            map <- getMapScore(wgdata, LR1)
            map <- aggregate(map, list(id=idx), mean, simplify=TRUE)
            binsGR[map[, 1]]\$map <- map[, 2]
            len <- coverage(LR)
            len <- sum(len)
            binsGR[names(len)]\$len <- len
            FragmentendsL <- trim(resize(FragmentendsL, 200, fix="end"))
            FragmentendsR <- trim(resize(FragmentendsR, 200, fix="start"))
            ends <- c(FragmentendsL, FragmentendsR)
            ol <- findOverlaps(ends, binsGR, type = "within", select = "all")
            LR <- GRanges(names(binsGR[subjectHits(ol)]), ranges(ends[queryHits(ol)]))
            LR <- reduce(LR)
            idx <- as.character(seqnames(LR))
            LR <- GRanges(sub("-.*\$", "", idx), ranges(LR), idx=idx)
            LR\$gc <- getGC(genome, LR, as.prob=FALSE)
            LR\$wid <- width(LR)
            gc <- rowsum(as.data.frame(mcols(LR[, c("gc", "wid")])), group=LR\$idx)
            gc_content <- gc[, 1]/gc[, 2]
            names(gc_content) <- rownames(gc)
            binsGR[names(gc_content)]\$gc <- gc_content
        }else{
            binsGR\$len <- NULL
            binsGR\$gc <- getGC(genome, binsGR)
            binsGR\$map <- getMapScore(wgdata, binsGR)
        }
        bintolen <- as.data.frame(binsGR)
        bintolen\$bins <- names(binsGR)
        bintolen <- bintolen[, c("bins", "gc", "map", "len"), drop=FALSE]
        if (sum(bintolen\$map) == 0) {
            bintolen <- bintolen %>% dplyr::select(-.data\$map)
        }
    }
    data.table::fwrite(bintolen, bintolenoutput, row.names = FALSE,
        quote = FALSE, sep = "\\t")
    """
}
