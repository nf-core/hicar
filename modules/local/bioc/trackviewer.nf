process BIOC_TRACKVIEWER {
    tag "$bin_size"
    label 'process_high'
    label 'process_long'
    label 'error_ignore'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(bin_size), path(cools), path(events), path(anchors)
    path raw_pairs // .unselected.pairs.gz of samfrag
    path gtf
    path chrom_sizes
    path restrict

    output:
    path "${prefix}/*"            , emit: v4c
    path "versions.yml"           , emit: versions

    script:
    prefix   = task.ext.prefix ?: "trackViewer_bin${bin_size}"
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 for trackViewer parser
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################

    pkgs <- c("trackViewer", "GenomicFeatures", "InteractionSet", "rtracklayer", "rhdf5")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    # Options
    ## make_option(c("-g", "--gtf"), type="character", default=NULL, help="filename of gtf", metavar="string")
    ## make_option(c("-s", "--chromsize"), type="character", default=NULL, help="filename of chrome size", metavar="string")
    ## make_option(c("-x", "--restrict"), type="character", default=NULL, help="filename of restrict cut", metavar="string")
    ## make_option(c("-r", "--resolution"), type="integer", default=NULL, help="resolution", metavar="integer")
    ## make_option(c("-d", "--gap"), type="integer", default=NULL, help="gap, default 2*resolution", metavar="integer")
    ## make_option(c("-l", "--readlength"), type="integer", default=NULL, help="reads length used to strengthen the signal, default resolution/10", metavar="integer")
    ## make_option(c("-m", "--maxevents"), type="integer", default=25, help="max events to plot", metavar="integer")
    ## make_option(c("-o", "--output"), type="character", default=".", help="output folder", metavar="string")
    ## make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer")
    ## make_option(c("-e", "--events"), type="character", default=NULL, help="given events csv file, must be ginteractions file", metavar="string")
    maxRegionWidth <- 1e6
    maxEvent <- 25
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
    option_list <- list("gap"=c("--gap", "-d", "integer"),
                        "readlength"=c("--readlength", "-l", "integer"),
                        "maxevents"=c("--maxevents", "-m", "integer"))
    opt <- parse_args(option_list, args)

    gtf <- "${gtf}"
    resolution <- ${bin_size}
    gap <- 2 * resolution
    readwidth <- ceiling(resolution/10)
    output <- "${prefix}"
    chrom_size <- "${chrom_sizes}"
    restrict_cut <- "${restrict}"
    if(!is.null(opt\$gap)){
        gap <- as.numeric(opt\$gap)
    }
    if(!is.null(opt\$readlength)){
        readwidth <- as.numeric(opt\$readlength)
    }
    if(!is.null(opt\$maxevents)){
        maxEvent <- as.numeric(opt\$maxevents)
    }

    # read the signals
    cools <- dir(".", "cool\$", full.names = TRUE, recursive = TRUE)
    pairfiles <- dir(".", ".h5\$", full.names = TRUE, recursive = TRUE)
    if(!is.null(opt\$events)){
        evts <- read.csv(opt\$events)
        if(!"fdr" %in% colnames(evts)){
            evts\$fdr <- rep(1, nrow(evts))
        }
    }else{
        evts <- dir(".", ".bedpe\$", full.names = TRUE, recursive = TRUE)
        if(length(evts)<1) stop("No events file.")
        header <- lapply(evts, read.delim, header=FALSE, nrow=1)
        evts <- mapply(evts, header, FUN=function(f, h){
            hasHeader <- all(c("chr1", "start1", "end1", "chr2", "start2", "end2") %in%
                                h[1, , drop=TRUE])
            .ele <- read.delim(f, header = hasHeader)
            if(!hasHeader){
                colnames(.ele)[1:6] <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
                ## bedpe, [start, end)
                .ele\$start1 <- .ele\$start1+1
                .ele\$start2 <- .ele\$start2+1
            }
            .ele
        }, SIMPLIFY = FALSE)
        evts <- do.call(rbind, evts)
    }
    stopifnot(nrow(evts)>0)
    colnames(evts) <- tolower(colnames(evts))
    stopifnot(all(c("chr1", "chr2", "start1", "start2", "end1", "end2") %in%
        colnames(evts)))
    if(!"fdr" %in% colnames(evts)){
        ## count the reads, TODO
        evts\$fdr <- 1
    }
    evts <- evts[order(evts\$fdr), , drop=FALSE]
    if(nrow(evts)<1) stop("No events in files.")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)

    ## get events ranges
    restrict_pos <- import(restrict_cut, format = "BED")
    seqlen <- read.delim(chrom_size, header=FALSE, row.names=1)
    seql <- as.numeric(seqlen[, 1])
    names(seql) <- rownames(seqlen)
    txdb <- makeTxDbFromGFF(gtf)
    gtf <- import(gtf, format = "GTF")
    map <- gtf\$gene_name
    names(map) <- gtf\$gene_id
    map <- map[!duplicated(names(map))]
    grs <- with(evts, GInteractions(anchor1 = GRanges(chr1, IRanges(start1, end1)),
                                    anchor2 = GRanges(chr2, IRanges(start2, end2)),
                                    mode = "strict",
                                    score = fdr))
    grs\$score <- 1-grs\$score
    grs <- unique(grs)
    grs_narrow <- grs
    regions(grs_narrow) <-
        resize(regions(grs), width = ceiling(width(regions(grs))/2), fix = "center")
    link <- gi2track(grs_narrow)
    setTrackStyleParam(link, "tracktype", "link")
    setTrackStyleParam(link, "color", c("gray80", "yellow", "brown"))
    reg0 <- GRanges(seqnames(first(grs)),
                            IRanges(start = start(first(grs)),
                                    end = end(second(grs))))
    reg0 <- reg0[width(reg0)<maxRegionWidth]
    reg <- reduce(reg0, min.gapwidth = gap)
    gr1 <- unique(subsetByOverlaps(reg, grs))
    gr1 <- gr1[seq.int(min(maxEvent, length(gr1)))]
    prettyMax <- function(x){
        if(x<1) return(round(x, digits = 2))
        if(x<=5) return(ceiling(x))
        if(x<10) return(2*ceiling(x/2))
        if(x<100) return(10*ceiling(x/10))
        if(x<1000) return(50*ceiling(x/50))
        n <- 10^ceiling(log10(x))
        return(n*ceiling(x/n))
    }

    ## read file and summary counts
    getIndex <- function(pos, tileWidth, ext=150){
        A <- ceiling((start(pos)-ext)/rep(tileWidth, length(pos)))
        B <- ceiling((end(pos)+ext)/rep(tileWidth, length(pos)))
        out <- mapply(A, B, FUN=seq, SIMPLIFY=FALSE)
        sort(unique(unlist(out)))
    }
    getPath <- function(root, ...){
        paste(root, ..., sep="/")
    }
    readPairs <- function(pair, chrom, range){
        tileWidth <- h5read(pair, "header/tile_width")
        idx <- getIndex(range, tileWidth)
        idx <- expand.grid(idx, idx)
        idx <- paste(idx[, 1], idx[, 2], sep="_")
        inf <- H5Fopen(pair, flags="H5F_ACC_RDONLY")
        on.exit(H5Fclose(inf))
        pc <- lapply(idx, function(.ele){
            n <- getPath("data", chrom, chrom, .ele, "position")
            if(H5Lexists(inf, n)){
                h5read(inf, n)
            }
        })
        pc <- do.call(rbind, pc)
        strand <- lapply(idx, function(.ele){
            n <- getPath("data", chrom, chrom, .ele, "strand")
            if(H5Lexists(inf, n)){
                h5read(inf, n)
            }
        })
        strand <- do.call(rbind, strand)
        H5Fclose(inf)
        h5closeAll()
        on.exit()
        GInteractions(anchor1 = GRanges(chrom, IRanges(pc[, 1], width=readwidth), strand = strand[, 1]),
                    anchor2 = GRanges(chrom, IRanges(pc[, 2], width=readwidth), strand = strand[, 2]))
    }
    loadPairFile <- function(filenames, ranges, resolution){
        stopifnot(is.character(filenames))
        stopifnot(all(file.exists(filenames)))
        start(ranges) <- start(ranges) - 2*resolution
        end(ranges) <- end(ranges) + 2*resolution
        names(filenames) <- sub(".h5\$", "", basename(filenames))
        total <- lapply(filenames, h5read, name="header/total")
        chrom <- as.character(seqnames(ranges)[1])
        ranges <- ranges[seqnames(ranges)==chrom]
        gi <- lapply(filenames, readPairs, chrom=chrom, range=ranges)
        list(gi=gi, total=total)
    }

    readPairFile <- function(filenames, ranges, resolution){
        stopifnot(is(ranges, "GRanges"))
        chunks <- loadPairFile(filenames, ranges, resolution)
        out <- list()
        total <- list()
        for(fn in names(chunks[["gi"]])){
            chunk <- chunks[["gi"]][[fn]]
            if(length(chunk)){
                total[[fn]] <- chunks[["total"]][[fn]]
                out[[fn]] <-
                    subsetByOverlaps(chunk, ranges = ranges, use.region="both")
            }
            rm(chunk)
        }
        ## split by group
        f <- sub("_REP\\\\d+", "", names(out))
        out <- split(out, f)
        giRbind <- function(a, b){
            GInteractions(anchor1 = c(first(a), first(b)),
                        anchor2 = c(second(a), second(b)),
                        regions=sort(unique(c(first(a), first(b), second(a), second(b)))))
        }
        out <- lapply(out, function(.ele){
            Reduce(giRbind, .ele)
        })
        total <- split(unlist(total), f)
        total <- sapply(total, sum)
        total <- median(total)/total
        return(list(gi=out, norm_factor=total))
    }

    for(i in seq_along(gr1)){
        if(i < maxEvent){
            try_res <- try({
                gr <- gr1[i]
                start(gr) <- start(gr) - 2* resolution
                end(gr) <- end(gr) + 2* resolution
                seqlevelsStyle(gr) <- seqlevelsStyle(txdb)[1]
                seqlengths(gr) <- seql[seqlevels(gr)]
                gr <- GenomicRanges::trim(gr)
                chr_gr <- as(seqinfo(gr), "GRanges")[seqnames(gr)[1]]
                bait <- subsetByOverlaps(grs, gr, use.region="both")
                ## merge the bait by ends
                bait <- split(second(bait), as.character(first(bait)))
                bait <- lapply(bait, range)
                bait <- unlist(GRangesList(bait))
                bait <- GInteractions(anchor1 = parse2GRanges(names(bait)),
                                    anchor2 = unname(bait))
                ## get the v4c
                info <- readPairFile(pairfiles,
                                    ranges = c(first(bait), second(bait)),
                                    resolution=resolution)
                h5closeAll()
                total <- info\$norm_factor
                names(cools) <- sub("\\\\d+\\\\.mcool\$", "", basename(cools))
                gis <- lapply(cools, importGInteractions,
                            resolution = resolution,
                            format="cool",
                            ranges=gr, out = "GInteractions")
                gis <- mapply(gis, total[names(gis)], FUN=function(.ele, .total){
                    .ele\$score <- log2(.ele\$score+1) * .total
                    .ele
                })
                maxV <- max(sapply(gis, function(.ele) max(.ele\$score, na.rm=TRUE)))
                maxV <- prettyMax(maxV)
                heat <- lapply(gis, function(.ele){
                    .ele <- gi2track(.ele)
                    setTrackStyleParam(.ele, "breaks", seq(from=0, to=maxV, by=maxV/10))
                    setTrackStyleParam(.ele, "color", c("lightblue", "yellow", "red"))
                    #setTrackStyleParam(.ele, "ylim", c(0, .5))
                    .ele
                })
                names(heat) <- paste("valid pairs", names(heat))
                get_v4c <- function(bait1, bait2){
                    v4c <- lapply(seq_along(bait1), function(.e){
                        vp1 <- mapply(info\$gi, total[names(info\$gi)], FUN=function(.ele, .total){
                            .ele <- subsetByOverlaps(.ele, bait1[.e], use.region="both")
                            .ele <- GRanges(coverage(c(first(.ele), second(.ele))))
                            .ele\$score <- .ele\$score*.total
                            new("track", dat=.ele,
                                type="data", format = "BED",
                                name = paste0(seqnames(bait1[.e]), ":",
                                            start(bait1[.e]), "-",
                                            end(bait1[.e])))
                        })
                        maxY <- sapply(vp1, FUN = function(.ele){
                            max(subsetByOverlaps(.ele\$dat, bait2[.e])\$score)
                        })
                        maxY <- prettyMax(max(maxY))
                        vp1 <- lapply(vp1, function(.ele){
                            setTrackStyleParam(.ele, "ylim", c(0, maxY))
                            setTrackStyleParam(.ele, "color", "gray80")
                            .ele
                        })
                    })
                    names(v4c) <- paste0(seqnames(bait1), ":",
                                        start(bait1), "-",
                                        end(bait1))
                    v4c <- unlist(v4c, recursive = FALSE)
                }
                v4c_first <- get_v4c(first(bait), second(bait))
                v4c_second <- get_v4c(second(bait), first(bait))
                ids <- getGeneIDsFromTxDb(gr, txdb)
                if(length(ids)){
                    genes <- geneTrack(ids, txdb, map[ids], asList=FALSE)
                    tL <- trackList(genes, link, heat, v4c_first, v4c_second,
                                    heightDist = c(1, 1, 2*length(heat),
                                    2*length(v4c_first), 2*length(v4c_second)))
                    names(tL)[2] <- "called links"
                }else{
                    tL <- trackList(link, heat, v4c_first, v4c_second,
                                    heightDist = c(1, 2*length(heat),
                                    2*length(v4c_first), 2*length(v4c_second)))
                    names(tL)[1] <- "called links"
                }
                pdf(file.path(output, paste0("event_", i, "_", seqnames(gr), ":", start(gr), "-", end(gr), ".pdf")),
                    width = 9, height = length(tL))
                viewTracks(tL, gr=gr, autoOptimizeStyle = TRUE)
                dev.off()
            })
            if(inherits(try_res, "try-error")){
                message(try_res)
            }
        }
    }
    """
}
