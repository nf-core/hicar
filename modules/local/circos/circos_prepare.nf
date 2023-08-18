process CIRCOS_PREPARE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(meta), path(bedpe)
    val ucscname
    path gtf
    path chromsize

    output:
    tuple val(meta), path("${prefix}/*")           , emit: circos
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Oct. 16, 2021 prepare data for circos
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    options(scipen=10)
    library(rtracklayer)
    versions <- c(
        "${task.process}:",
        paste("    rtracklayer:", as.character(packageVersion("rtracklayer"))))
    writeLines(versions, "versions.yml")

    ## Options
    chromsize <- "$chromsize"
    gtf <- "$gtf"
    ucscname <- "$ucscname"
    outfolder <- "${prefix}"
    totalLinks <- 1e4

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
    option_list <- list("totalLinks"=c("--totalLinks", "-n", "numeric"))
    opt <- parse_args(option_list, args)
    if(!is.null(opt[["totalLinks"]])){
        totalLinks <- opt[["totalLinks"]]
    }

    getScore <- function(seql, rg){
        gtile <- tileGenome(seqlengths = seql, tilewidth = 1e5)
        gtile <- unlist(gtile)
        gtile <- split(gtile, seqnames(gtile))
        sharedChr <- intersect(names(rg), names(gtile))
        rg <- rg[sharedChr]
        gtile <- gtile[names(rg)]
        vw <- Views(rg, gtile)
        vm <- viewSums(vw, na.rm = TRUE)
        stopifnot(identical(names(vm), names(gtile)))
        gtile <- unlist(gtile, use.names = FALSE)
        gtile\$score <- unlist(vm, use.names = FALSE)
        gtile
    }
    filterSeqnames <- function(x){
        seqlevelsStyle(x) <- "UCSC"
        x <- x[seqnames(x) %in% standardChromosomes(x)]
        x <- keepStandardChromosomes(x)
        x
    }
    dir.create(outfolder, showWarnings = FALSE)
    gtf <- import(gtf)
    ## create karyotype file
    chromsize <- read.delim(chromsize, header=FALSE)
    chromsize[, 1] <- as.character(chromsize[, 1])
    chromsize <- chromsize[chromsize[, 1] %in% standardChromosomes(gtf), , drop=FALSE]
    seqlevelsStyle(chromsize[, 1]) <- "UCSC"
    emptyTADfile <- c(paste(chromsize[, 1], 1, 1, 0), paste(chromsize[, 1], 2, 2, 1))
    chromsize <- cbind("chr", "-", chromsize[, c(1, 1)], 0, chromsize[, 2],
                        paste0("chr", sub("^chr", "", chromsize[, 1])))
    colnames(chromsize) <- c("chr", "-", "chrname", "chrlabel",
                            0, "chrlen", "chrcolor")
    chromsize <- chromsize[order(as.numeric(sub("chr", "", chromsize[, "chrname"])),
                            chromsize[, "chrname"]), , drop=FALSE]
    write.table(chromsize, file.path(outfolder, "karyotype.tab"),
                quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")

    ## create label file
    seqn <- sort(as.character(unique(chromsize\$chrname)))[1]
    seql <- chromsize\$chrlen
    names(seql) <- chromsize\$chrname
    labelA <- labelB <- c(seqn, 0, 5000)
    labelA <- c(labelA, "TAD")
    labelB <- c(labelB, "compartments")
    writeLines(labelA, file.path(outfolder, "labelA.txt"), sep=" ")
    writeLines(labelB, file.path(outfolder, "labelB.txt"), sep=" ")

    (interaction <- dir(".", pattern = ".bedpe", ignore.case = TRUE))
    (compartment <- dir(".", pattern = ".bigWig|bw\$", ignore.case = TRUE))
    (tad <- dir(".", pattern = ".bed\$", ignore.case = TRUE))
    writeLines(character(0L), file.path(outfolder, "link.txt"))
    try({
        headerline <- readLines(interaction[1], n=1)
        if(grepl("start1", headerline[1])){
            ## output from MAPS
            pe <- read.delim(interaction)
            pe <- Pairs(GRanges(pe[, "chr1"], IRanges(pe[, "start1"]+1, pe[, "end1"])),
                        GRanges(pe[, "chr2"], IRanges(pe[, "start2"]+1, pe[, "end2"])),
                        score = pe[, "ClusterNegLog10P"])
        }else{
            pe <- import(interaction, format="BEDPE")
        }
        seqlevelsStyle(first(pe)) <- seqlevelsStyle(second(pe)) <- "UCSC"
        keep <- seqnames(first(pe)) %in% standardChromosomes(first(pe)) &
            seqnames(second(pe)) %in% standardChromosomes(second(pe))
        pe <- pe[keep]
        pes <- unique(pe[order(mcols(pe)\$score, decreasing=TRUE)])
        if(length(pes)>totalLinks){
            ## keep top 24K links only. otherwise hard to plot.
            pes <- pes[seq.int(min(totalLinks, length(pes)))]
        }
        out <- as.data.frame(pes)
        scores <- sqrt(range(mcols(pe)\$score)/10)
        scores <- c(floor(scores[1]), ceiling(scores[2]))
        cid <- cut(sqrt(mcols(pes)\$score/10), breaks = seq(scores[1], scores[2]))
        levels(cid) <- seq_along(levels(cid))
        out <- cbind(out[, c("first.seqnames", "first.start", "first.end",
                            "second.seqnames", "second.start", "second.end")],
                    thickness=paste0("thickness=", cid))

        write.table(out, file.path(outfolder, "link.txt"),
            quote=FALSE, col.names=FALSE, row.names=FALSE,
            sep=" ") ## output cis- and trans- interactions
    })

    ##tads
    writeLines(emptyTADfile, file.path(outfolder, "tad.txt"))
    try({
        rg <- tryCatch({
            rg <- import(tad[1], format="BED")
        }, error=function(.e){
            rg <- read.table(tad[1])
            rg <- rg[, c(seq.int(3), 8)]
            colnames(rg) <- c("seqnames", "start", "end", "score")
            rg <- GRanges(rg)
        })
        rg <- filterSeqnames(rg)
        rg <- sort(rg)
        rg <- as.data.frame(rg)
        rg <- rg[, c("seqnames", "start", "end", "score"), drop=FALSE]
        write.table(rg, file.path(outfolder, "tad.txt"),
                    quote=FALSE, col.names=FALSE, row.names=FALSE,
                    sep=" ")
    })
    ## Compartment
    writeLines(character(0L), file.path(outfolder, "compartment.txt"))
    try({
        cp <- import(compartment[1], format="bigWig")
        cp <- filterSeqnames(cp)
        cp\$score[is.na(cp\$score)] <- 0
        if(all(cp\$score==0)){
            cp\$score <- rep(c(-1, 1), length(cp))[seq_along(cp)]
        }
        cp\$score <- ifelse(cp\$score==0, 0, ifelse(cp\$score<0, -1, 1))
        if(sum(cp\$score==-1)==0){
            cp\$score[cp\$score==0] <- -1
        }
        if(sum(cp\$score==-1)==0){
            gaps <- gaps(cp)
            gaps <- gaps[strand(gaps)=="*"]
            gaps\$score <- -1
            cp <- sort(c(cp, gaps))
        }
        cvg <- coverage(cp, weight=cp\$score)
        cp <- as(cvg, "GRanges")
        cp <- cp[cp\$score!=0]
        cp <- as.data.frame(cp)
        cp <- cp[, c("seqnames", "start", "end", "score"), drop=FALSE]
        write.table(cp, file.path(outfolder, "compartment.txt"),
                    quote=FALSE, col.names = FALSE, row.names = FALSE,
                    sep=" ")
    })

    ## create karyotype file again
    try_res <- try({
        session <- browserSession()
        genome(session) <- ucscname
        ideo <- getTable(ucscTableQuery(session, table="cytoBandIdeo"))
        ideo <- ideo[ideo\$chrom %in% chromsize\$chrname, , drop=FALSE]
        ideo <- data.frame("band", ideo\$chrom, ideo\$name, ideo\$name,
                    ideo\$chromStart, ideo\$chromEnd, ideo\$gieStain)
        colnames(ideo) <- colnames(chromsize)
        ideo <- rbind(chromsize, ideo)
        for(i in c("chrname", "chrlabel")){
            ideo[ideo[, i]=="", i] <- make.names(ideo[ideo[, i]=="", 2], unique=TRUE)
        }
        write.table(ideo, file.path(outfolder, "karyotype.tab"),
                    quote=FALSE, col.names=FALSE, row.names=FALSE, sep=" ")
    })
    if(inherits(try_res, "try-error")){
        message(try_res)
    }
    """
}
