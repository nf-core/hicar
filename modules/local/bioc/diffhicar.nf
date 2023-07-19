process DIFFHICAR {
    tag "$bin_size"
    label 'process_medium'

    conda "bioconda::bioconductor-edger=3.36.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-edger:3.36.0--r41hc247a5b_2' :
        'biocontainers/bioconductor-edger:3.36.0--r41hc247a5b_2' }"

    input:
    tuple val(bin_size), path(peaks, stageAs: "peaks/*"), path(long_bedpe, stageAs: "long/*")
    val long_bedpe_postfix

    output:
    tuple val(bin_size), path("${prefix}/*")               , emit: diff
    tuple val(bin_size), val("$prefix"), path("${prefix}/edgeR.DEtable*") , optional: true, emit: anno
    tuple val(bin_size), val("$prefix"), path("${prefix}/*.bedpe") , optional: true, emit: bedpe
    path "${prefix}/*.qc.json"                             , emit: stats
    path "versions.yml"                                    , emit: versions

    script:
    prefix   = task.ext.prefix ?: "edgeR_bin${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 call edgeR
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    library(edgeR)
    versions <- c(
        "${task.process}:",
        paste("    edgeR:", as.character(packageVersion("edgeR"))))
    writeLines(versions, "versions.yml")

    binsize = "$prefix"

    ## get peaks
    pf <- dir("peaks", "bedpe", full.names = TRUE)
    header <- lapply(pf, read.delim, header=FALSE, nrow=1)
    peaks <- mapply(pf, header, FUN=function(f, h){
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
    ### reduce the peaks
    peaks <- unique(do.call(rbind, peaks)[, c("chr1", "start1", "end1",
                                            "chr2", "start2", "end2")])

    ## get counts
    pc <- dir("long", "bedpe", full.names = FALSE)
    cnts <- lapply(file.path("long", pc), read.table)
    samples <- sub("(_REP\\\\d+)\\\\.(.*?)\\\\.${long_bedpe_postfix}", "\\\\1", pc)
    cnts <- lapply(split(cnts, samples), do.call, what=rbind)
    sizeFactor <- vapply(cnts, FUN=function(.ele) sum(.ele[, 7], na.rm = TRUE),
                        FUN.VALUE = numeric(1))

    ## check coordinates style
    ### bedpe, [start, end)
    checkStyle <- function(x, y){
        if(sum((x[, 2]+1) %in% y[, 2]) > sum(x[, 2] %in% y[, 2])) return(1)
        if(sum((y[, 2]+1) %in% x[, 2]) > sum(y[, 2] %in% x[, 2])) return(-1)
        return(0)
    }
    offset <- checkStyle(peaks, do.call(rbind, cnts))
    peaks\$start1 <- peaks\$start1 + offset
    peaks\$start2 <- peaks\$start2 + offset
    getID <- function(mat) gsub("\\\\s+", "", apply(mat[, seq.int(6)], 1, paste, collapse="_"))
    getID1 <- function(mat) gsub("\\\\s+", "", apply(mat[, seq.int(3)], 1, paste, collapse="_"))
    getID2 <- function(mat) gsub("\\\\s+", "", apply(mat[, 4:6], 1, paste, collapse="_"))
    ## prefilter, to decrease the memory cost
    peaks_id1 <- getID1(peaks)
    peaks_id2 <- getID2(peaks)
    cnts <- lapply(cnts, function(.ele) .ele[getID1(.ele) %in% peaks_id1, , drop=FALSE])
    cnts <- lapply(cnts, function(.ele) .ele[getID2(.ele) %in% peaks_id2, , drop=FALSE])
    rm(peaks_id1, peaks_id2, getID1, getID2)
    ## match all the counts for peaks
    peaks_id <- getID(peaks)
    cnts <- do.call(cbind, lapply(cnts, function(.ele){
        .ele[match(peaks_id, getID(.ele)), 7]
    }))
    cnts[is.na(cnts)] <- 0
    names(peaks_id) <- paste0(rep("p", length(peaks_id)), seq_along(peaks_id))
    rownames(cnts) <- names(peaks_id)
    if(all(colSums(cnts)==0)){
        stop("Can not get counts table.")
    }

    pf <- as.character(binsize)
    dir.create(pf, showWarnings = FALSE, recursive=TRUE)

    fname <- function(subf, ext, ...){ # create file name
        pff <- ifelse(is.na(subf), pf, file.path(pf, subf))
        dir.create(pff, showWarnings = FALSE, recursive = TRUE)
        file.path(pff, paste(..., ext, sep="."))
    }

    ## write counts
    write.csv(cbind(peaks, cnts), fname(NA, "csv", "raw.counts"), row.names = FALSE)
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

    if(length(contrasts.lev)>1 && any(table(condition)>1)){
        contrasts <- combn(contrasts.lev, 2, simplify = FALSE) ## pair all conditions
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
        try_res <- try({ ## try to output PCA results for multiQC
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
            writeLines(json, fname(NA, "json", "Multidimensional.scaling.qc"))
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
            res <- cbind(peaks, res[names(peaks_id), ])
            write.csv(res, fname(name, "csv", "edgeR.DEtable", name), row.names = FALSE)
            ## save metadata
            elementMetadata <- do.call(rbind, lapply(c("adjust.method","comparison","test"), function(.ele) rs[[.ele]]))
            rownames(elementMetadata) <- c("adjust.method","comparison","test")
            colnames(elementMetadata)[1] <- "value"
            write.csv(elementMetadata, fname(name, "csv", "edgeR.metadata", name), row.names = TRUE)
            ## save subset results
            res.s <- res[res\$FDR<0.05 & abs(res\$logFC)>1, ]
            write.csv(res.s, fname(name, "csv", "edgeR.DEtable", name, "padj0.05.lfc1"), row.names = FALSE)
            ## Volcano plot
            res\$qvalue <- -10*log10(res\$PValue)
            res.s\$qvalue <- -10*log10(res.s\$PValue)
            res.s1 <- cbind(res.s[, c("chr1", "start1", "end1", "chr2", "start2", "end2", "qvalue"), drop = FALSE],
                name = rep('.', nrow(res.s)),
                res.s[, c("qvalue"), drop = FALSE],
                strand1 = rep('.', nrow(res.s)), strand2 = rep('.', nrow(res.s)))
            write.table(res.s1, fname(name, "bedpe", "edgeR.DEtable", name, "padj0.05.lfc1"), row.names=FALSE, col.names=FALSE, sep="\\t", quote=FALSE)
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
