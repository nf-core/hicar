process DIFFHIC {
    tag "$bin_size"
    label 'process_medium'

    conda "bioconda::bioconductor-diffhic=1.26.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-diffhic:1.26.0--r41hc247a5b_2' :
        'biocontainers/bioconductor-diffhic:1.26.0--r41hc247a5b_2' }"

    input:
    tuple val(bin_size), path(peaks, stageAs: "peaks/*"), path(long_bedpe, stageAs: "long/*")
    val long_bedpe_postfix

    output:
    tuple val(bin_size), path("${prefix}/*")                                               , emit: diff
    tuple val(bin_size), val("$prefix"), path("${prefix}/diffHic.DEtable*"), optional: true, emit: anno
    tuple val(bin_size), val("$prefix"), path("${prefix}/*.bedpe")         , optional: true, emit: bedpe
    path "${prefix}/*.qc.json"                                                             , emit: stats
    path "versions.yml"                                                                    , emit: versions

    script:
    prefix   = task.ext.prefix ?: "diffHic_bin${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Nov. 11, 2022 call diffhic
    ## Copyright (c) 2022 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("edgeR", "diffHic")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    OUTFOLDER = "$prefix"

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
    peaks <- with(peaks, GInteractions(GRanges(chr1, IRanges(start1, end1)),
        GRanges(chr2, IRanges(start2, end2))))
    ## get counts
    pc <- dir("long", "bedpe", full.names = FALSE)
    cnts <- lapply(file.path("long", pc), read.table)
    samples <- sub("(_REP\\\\d+)\\\\.(.*?)\\\\.${long_bedpe_postfix}", "\\\\1", pc)
    cnts <- lapply(split(cnts, samples), do.call, what=rbind)
    sizeFactor <- vapply(cnts, FUN=function(.ele) sum(.ele[, 7], na.rm = TRUE),
        FUN.VALUE = numeric(1))
    cnts <- lapply(cnts, function(.ele){
    with(.ele, GInteractions(GRanges(V1, IRanges(V2, V3)),
        GRanges(V4, IRanges(V5, V6)),
        score = V7))
    })
    cnts0 <- cnts ## backup
    cnts <- lapply(cnts, function(.ele){
        ol <- findOverlaps(peaks, .ele, type="equal", select = "first")
        .peak <- peaks
        .peak\$score <- 0
        .peak\$score[!is.na(ol)] <- .ele[ol[!is.na(ol)]]\$score
        .peak\$score
    })
    cnts <- do.call(cbind, cnts)
    if(all(colSums(cnts)==0)){## peaks not exactly bin based
        cnts <- lapply(cnts0, function(.ele){
            ol <- findOverlaps(peaks, .ele, select = "all")
            .peak <- peaks
            .peak\$score <- 0
            .score <- lapply(split(.ele[subjectHits(ol)]\$score, queryHits(ol)), max, na.rm=TRUE)
            .score <- unlist(.score)
            .peak\$score[as.numeric(names(.score))] <- .score
            .peak\$score
        })
        cnts <- do.call(cbind, cnts)
    }

    pf <- as.character(OUTFOLDER)
    dir.create(pf, showWarnings = FALSE, recursive=TRUE)

    fname <- function(subf, ext, ...){ # create file name
        pff <- ifelse(is.na(subf), pf, file.path(pf, subf))
        dir.create(pff, showWarnings = FALSE, recursive = TRUE)
        file.path(pff, paste(..., ext, sep="."))
    }

    ## write counts
    write.csv(cbind(as.data.frame(peaks), cnts), fname(NA, "csv", "raw.counts"), row.names = FALSE)
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

    data <- InteractionSet(assays = list(counts=cnts),
        interactions = peaks,
        colData = DataFrame(totals = sizeFactor))
    data.tmm <- csaw::normOffsets(data, se.out=TRUE)
    pdf(fname(NA, "pdf", "normalization"), width=12, height=6)
    par(mfrow=c(1, 2))
    ab <- aveLogCPM(asDGEList(data))
    o <- order(ab)
    adj.counts <- cpm(asDGEList(data), log=TRUE)
    mval <- adj.counts[,3]-adj.counts[,2]
    plotMD(asDGEList(data), main="before")
    fit <- loessFit(x=ab, y=mval)
    lines(ab[o], fit\$fitted[o], col="red")
    abline(h=0, col='blue', lty=2)

    ab <- aveLogCPM(asDGEList(data.tmm))
    o <- order(ab)
    adj.counts <- cpm(asDGEList(data.tmm), log=TRUE)
    mval <- adj.counts[,3]-adj.counts[,2]
    plotMD(asDGEList(data.tmm), main="after")
    fit <- loessFit(x=ab, y=mval)
    lines(ab[o], fit\$fitted[o], col="red")
    abline(h=0, col='blue', lty=2)

    dev.off()

    if(length(contrasts.lev)>1 && any(table(condition)>1)){
        contrasts <- combn(contrasts.lev, 2, simplify = FALSE) ## pair all conditions
        ## create DGEList
        group <- coldata\$condition
        y <- asDGEList(data.tmm,
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
            res <- cbind(peaks, res)
            write.csv(res, fname(name, "csv", "diffHic.DEtable", name), row.names = FALSE)
            ## save metadata
            elementMetadata <- do.call(rbind, lapply(c("adjust.method","comparison","test"), function(.ele) rs[[.ele]]))
            rownames(elementMetadata) <- c("adjust.method","comparison","test")
            colnames(elementMetadata)[1] <- "value"
            write.csv(elementMetadata, fname(name, "csv", "diffHic.metadata", name), row.names = TRUE)
            current_peaks <- peaks
            current_peaks\$score <- res\$F
            current_peaks <- current_peaks[res\$FDR<0.05 & abs(res\$logFC)>1]
            rtracklayer::export(current_peaks, fname(name, "bedpe", "diffHic.DEtable", name, "padj0.05.lfc1"), format = 'bedpe')
            ## save subset results
            res.s <- res[res\$FDR<0.05 & abs(res\$logFC)>1, ]
            write.csv(res.s, fname(name, "csv", "diffHic.DEtable", name, "padj0.05.lfc1"), row.names = FALSE)
            ## Volcano plot
            res\$qvalue <- -10*log10(res\$PValue)
            res.s\$qvalue <- -10*log10(res.s\$PValue)
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
