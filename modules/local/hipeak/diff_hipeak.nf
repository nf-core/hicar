process DIFF_HIPEAK {
    label 'process_medium'

    conda "bioconda::bioconductor-diffhic=1.24.0"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-diffhic:1.24.0--r41h399db7b_0 "
    } else {
        container "biocontainers/bioconductor-diffhic:1.24.0--r41h399db7b_0"
    }

    input:
    path peaks, stageAs: "peaks/*"
    path distalpair, stageAs: "pairs/*"

    output:
    path "${prefix}/*"                        , emit: diff
    path "${prefix}/*.bedpe"  , optional: true, emit: bedpe
    path "${prefix}/*.qc.json", optional: true, emit: stats
    path "versions.yml"                       , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix   = task.ext.prefix ? "${task.ext.prefix}" : "diffhicar"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 call edgeR
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("edgeR", "InteractionSet", "rhdf5", "BiocParallel")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

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
    option_list <- list("snow_type"=c("--snow_type", "-t", "character"))
    opt <- parse_args(option_list, strsplit("$args", "\\\\s+")[[1]])
    PREFIX <- "$prefix"
    NCORE <- as.numeric("$task.cpus")
    SNOW_TYPE <- "SOCK"
    if(!is.null(opt\$snow_type)){
        SNOW_TYPE <- opt\$snow_type
    }

    ## get peaks
    pf <- dir("peaks", "peaks", full.names = TRUE)
    peaks <- lapply(pf, read.table, header=TRUE)
    ### reduce the peaks
    peaks <- unique(do.call(rbind, peaks)[, c("chr1", "start1", "end1",
                                            "chr2", "start2", "end2")])
    peaks <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)),
                                        GRanges(chr2, IRanges(start2, end2))))
    reducePeaks <- function(x){
        y <- reduce(x)
        ol <- findOverlaps(x, y)
        stopifnot(all(seq_along(x) %in% queryHits(ol)))
        ol <- as.data.frame(ol)
        y[ol[match(seq_along(x), ol\$queryHits), "subjectHits"]]
    }
    first <- reducePeaks(first(peaks))
    second <- reducePeaks(second(peaks))
    peaks <- unique(GInteractions(first, second))

    ## get counts
    if(SNOW_TYPE=="FORK"){
        param <- MulticoreParam(workers = NCORE, progressbar = TRUE)
    }else{
        param <- SnowParam(workers = NCORE, progressbar = TRUE, type = SNOW_TYPE)
    }
    pc <- dir("pairs", "h5\$", full.names = FALSE)
    countByOverlaps <- function(pairs, peaks, sep="___"){
        getPath <- function(root, ...){
            paste(root, ..., sep="/")
        }
        readPairs <- function(pair, chrom1, chrom2){
            h5content <- rhdf5::h5ls(pair)
            h5content <- h5content[, "group"]
            h5content <- h5content[grepl("data.*\\\\d+_\\\\d+", h5content)]
            h5content <- unique(h5content)
            n <- h5content[grepl(paste0("data.", chrom1, ".", chrom2), h5content)]
            n <- getPath(n, "position")
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
        cnt <- lapply(names(peaks), function(chr){
            .peak <- peaks[[chr]]
            chr_ <- strsplit(chr, sep)[[1]]
            chrom1 <- chr_[1]
            chrom2 <- chr_[2]
            ps <- readPairs(pairs, chrom1, chrom2)
            counts_total <- rhdf5::h5read(pairs, "header/total")
            if(length(ps)<1){
                return(NULL)
            }
            ps <- InteractionSet::GInteractions(
                    GenomicRanges::GRanges(chrom1, IRanges::IRanges(ps[, 1], width=150)),
                    GenomicRanges::GRanges(chrom2, IRanges::IRanges(ps[, 2], width=150)))
            counts_tab <- IRanges::countOverlaps(.peak, ps, use.region="both")
            counts_tab <- cbind(ID=.peak\$ID, counts_tab)
            list(count=counts_tab, total=counts_total)
        })
        cnt <- cnt[lengths(cnt)>0]
        counts_total <- vapply(cnt, FUN=function(.ele) .ele\$total,
                            FUN.VALUE = numeric(1))
        counts_total <- sum(counts_total)
        counts_tab <- do.call(rbind, lapply(cnt, function(.ele) .ele\$count))
        list(count=counts_tab, total=counts_total)
    }

    peaks\$ID <- seq_along(peaks)
    peaks.s <- split(peaks, paste(seqnames(first(peaks)), seqnames(second(peaks)), sep="___"))
    try_res <- try({cnts <- bplapply(file.path("pairs", pc), countByOverlaps, peaks=peaks.s, sep="___", BPPARAM = param)})
    sizeFactor <- vapply(cnts, FUN=function(.ele) .ele\$total,
                        FUN.VALUE = numeric(1))
    if(inherits(try_res, "try-error") || all(sizeFactor==0)){ # check sizeFactor to make sure bplapply work
        cnts <- lapply(file.path("pairs", pc), countByOverlaps, peaks=peaks.s, sep="___")
    }
    h5closeAll()
    rm(peaks.s)
    samples <- sub("(_REP\\\\d+)\\\\.(.*?)h5\$", "\\\\1", pc)
    sizeFactor <- vapply(cnts, FUN=function(.ele) .ele\$total,
                        FUN.VALUE = numeric(1))
    names(sizeFactor) <- samples
    cnts <- lapply(cnts, function(.ele) .ele\$count)
    cnts <- mapply(cnts, samples, FUN=function(.d, .n){
        colnames(.d)[colnames(.d)!="ID"] <- .n
        .d
    }, SIMPLIFY=FALSE)
    cnts <- Reduce(function(x, y) merge(x, y, by="ID"), cnts)
    cnts <- cnts[match(peaks\$ID, cnts[, "ID"]), , drop=FALSE]
    cnts <- cnts[, colnames(cnts)!="ID", drop=FALSE]
    colnames(cnts) <- samples
    rownames(cnts) <- seq_along(peaks)
    mcols(peaks) <- cnts

    pf <- as.character(PREFIX)
    dir.create(pf)

    fname <- function(subf, ext, ...){
        pff <- ifelse(is.na(subf), pf, file.path(pf, subf))
        dir.create(pff, showWarnings = FALSE, recursive = TRUE)
        file.path(pff, paste(..., ext, sep="."))
    }

    ## write counts
    write.csv(peaks, fname(NA, "csv", "raw.counts"), row.names = FALSE)
    ## write sizeFactors
    write.csv(sizeFactor, fname(NA, "csv", "library.size"), row.names = TRUE)

    ## coldata
    sampleNames <- colnames(cnts)
    condition <- make.names(sub("_REP.*\$", "", sampleNames), allow_=TRUE)
    coldata <- data.frame(condition=factor(condition),
                        row.names = sampleNames)
    ## write designtable
    write.csv(coldata, fname(NA, "csv", "designTab"), row.names = TRUE)

    contrasts.lev <- levels(coldata\$condition)

    if(length(unique(contrasts.lev))>1 && any(table(condition)>1)){
        contrasts <- combn(contrasts.lev, 2, simplify = FALSE)
        ## create DGEList
        group <- coldata\$condition
        y <- DGEList(counts = cnts,
                    lib.size = sizeFactor,
                    group = group)

        ## do differential analysis
        names(contrasts) <- vapply(contrasts,
                                    FUN=paste,
                                    FUN.VALUE = character(1),
                                    collapse = "-")
        y <- calcNormFactors(y)
        design <- model.matrix(~0+group)
        colnames(design) <- levels(y\$samples\$group)
        y <- estimateDisp(y,design)
        fit <- glmQLFit(y, design)

        ## PCA
        pdf(fname(NA, "pdf", "Multidimensional.scaling.plot-plot"))
        mds <- plotMDS(y)
        dev.off()
        ## PCA for multiQC
        try_res <- try({
        json <- data.frame(x=mds\$x, y=mds\$y)
        rownames(json) <- colnames(y)
        json <- split(json, coldata[rownames(json), "condition"])
        json <- mapply(json, rainbow(n=length(json)), FUN=function(.ele, .color){
            .ele <- cbind(.ele, "name"=rownames(.ele))
            .ele <- apply(.ele, 1, function(.e){
                x <- names(.e)
                y <- .e
                .e <- paste0('{"x":', .e[1],
                            ', "y":', .e[2],
                            ', "color":"', .color,
                            '", "name":"', .e[3],
                            '"}')
            })
            .ele <- paste(.ele, collapse=", ")
        })
        json <- c(
                "{",
                '"id":"sample_pca",',
                '"data":{',
                paste(unlist(json), collapse=", "),
                "}",
                "}")
        writeLines(json, fname(NA, "json", "HiPeak.Multidimensional.scaling.qc"))
        })
        if(inherits(try_res, "try-error")){
            message(try_res)
        }

        ## plot dispersion
        pdf(fname(NA, "pdf", "DispersionEstimate-plot"))
        plotBCV(y)
        dev.off()
        ## plot QL dispersions
        pdf(fname(NA, "pdf", "Quasi-Likelihood-DispersionEstimate-plot"))
        plotQLDisp(fit)
        dev.off()

        res <- mapply(contrasts, names(contrasts), FUN = function(cont, name){
            BvsA <- makeContrasts(contrasts = name, levels = design)
            qlf <- glmQLFTest(fit, contrast = BvsA)
            rs <- topTags(qlf, n = nrow(qlf), sort.by = "none")
            ## MD-plot
            pdf(fname(name, "pdf", "Mean-Difference-plot", name))
            plotMD(qlf)
            abline(h=0, col="red", lty=2, lwd=2)
            dev.off()
            ## PValue distribution
            pdf(fname(name, "pdf", "PValue-distribution-plot", name))
            hist(rs\$table\$PValue, breaks = 20)
            dev.off()
            ## save res
            res <- as.data.frame(rs)
            res <- cbind(peaks[as.numeric(rownames(res))], res)
            colnames(res) <- sub("seqnames", "chr", colnames(res))
            write.csv(res, fname(name, "csv", "edgeR.DEtable", name), row.names = FALSE)
            ## save metadata
            elementMetadata <- do.call(rbind, lapply(c("adjust.method","comparison","test"), function(.ele) rs[[.ele]]))
            rownames(elementMetadata) <- c("adjust.method","comparison","test")
            colnames(elementMetadata)[1] <- "value"
            write.csv(elementMetadata, fname(name, "csv", "edgeR.metadata", name), row.names = TRUE)
            ## save subset results
            res.s <- res[res\$FDR<0.05 & abs(res\$logFC)>1, ]
            if(nrow(res.s)>0){
                write.csv(res.s, fname(name, "csv", "edgeR.DEtable", name, "padj0.05.lfc1"), row.names = FALSE)
                write.table(res.s[, c("chr1", "start1", "end1", "chr2", "start2", "end2"), drop=FALSE],
                    fname(name, "bedpe", "edgeR.DEtable", name, "padj0.05.lfc1"),
                    row.names = FALSE, col.names=FALSE, quote=FALSE, sep="\\t")
            }
            ## Volcano plot
            res\$qvalue <- -10*log10(res\$PValue)
            pdf(fname(name, "pdf", "Volcano-plot", name))
            plot(x=res\$logFC, y=res\$qvalue,
                main = paste("Volcano plot for", name),
                xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
                type = "p", col=NA)
            res.1 <- res
            if(nrow(res.1)>0) points(x=res.1\$logFC, y=res.1\$qvalue, pch = 20, cex=.5, col="gray80")
            if(nrow(res.s)>0) points(x=res.s\$logFC, y=res.s\$qvalue, pch = 19, cex=.5, col=ifelse(res.s\$logFC>0, "brown", "darkblue"))
            dev.off()
            res\$qvalue <- -10*log10(res\$PValue)
            png(fname(name, "png", "Volcano-plot", name))
            plot(x=res\$logFC, y=res\$qvalue,
                main = paste("Volcano plot for", name),
                xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
                type = "p", col=NA)
            res.1 <- res
            if(nrow(res.1)>0) points(x=res.1\$logFC, y=res.1\$qvalue, pch = 20, cex=.5, col="gray80")
            if(nrow(res.s)>0) points(x=res.s\$logFC, y=res.s\$qvalue, pch = 19, cex=.5, col=ifelse(res.s\$logFC>0, "brown", "darkblue"))
            dev.off()
        })
    }
    """
}
