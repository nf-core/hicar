#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on April. 29, 2021 call edgeR
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(edgeR)
writeLines(as.character(packageVersion("edgeR")), "edgeR.version.txt")

binsize = "diffhic_bin5000"
## load n.cores
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  binsize <- args[[1]]
  args <- args[-1]
}
if(length(args)>0){
  args <- lapply(args, function(.ele) eval(parse(.ele)))
}else{
  args <- list()
}


## get peaks
pf <- dir("peaks", "bedpe", full.names = TRUE)
peaks <- lapply(pf, read.delim)
### reduce the peaks
peaks <- unique(do.call(rbind, peaks)[, c("chr1", "start1", "end1",
                                          "chr2", "start2", "end2")])

## get counts
pc <- dir("long", "bedpe", full.names = FALSE)
cnts <- lapply(file.path("long", pc), read.table)
samples <- sub("(_REP\\d+)\\.(.*?)\\.long.intra.bedpe", "\\1", pc)
cnts <- lapply(split(cnts, samples), do.call, what=rbind)
sizeFactor <- vapply(cnts, FUN=function(.ele) sum(.ele[, 7], na.rm = TRUE),
                     FUN.VALUE = numeric(1))

getID <- function(mat) gsub("\\s+", "", apply(mat[, seq.int(6)], 1, paste, collapse="_"))
peaks_id <- getID(peaks)
cnts <- do.call(cbind, lapply(cnts, function(.ele){
  .ele[match(peaks_id, getID(.ele)), 7]
}))
cnts[is.na(cnts)] <- 0
names(peaks_id) <- paste0("p", seq_along(peaks_id))
rownames(cnts) <- names(peaks_id)

pf <- as.character(binsize)
dir.create(pf)

fname <- function(subf, ext, ...){
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
condition <- sub("_REP.*$", "", sampleNames)
coldata <- data.frame(condition=factor(condition),
                      row.names = sampleNames)
## write designtable
write.csv(coldata, fname(NA, "csv", "designTab"), row.names = TRUE)

contrasts.lev <- levels(coldata$condition)

if(length(contrasts.lev)>1 || any(table(condition)>1)){
  contrasts <- combn(contrasts.lev, 2, simplify = FALSE)
  ## create DGEList
  group <- coldata$condition
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
  colnames(design) <- levels(y$samples$group)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y, design)
  
  ## PCA
  pdf(fname(NA, "pdf", "Multidimensional.scaling.plot-plot"))
  plotMDS(y)
  dev.off()
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
    hist(rs$table$PValue, breaks = 20)
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
    res.s <- res[res$FDR<0.05 & abs(res$logFC)>1, ]
    write.csv(res.s, fname(name, "csv", "edgeR.DEtable", name, "padj0.05.lfc1"), row.names = FALSE)
    ## Volcano plot
    res$qvalue <- -10*log10(res$PValue)
    pdf(fname(name, "pdf", "Volcano-plot", name))
    plot(x=res$logFC, y=res$qvalue,
         main = paste("Volcano plot for", name),
         xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
         type = "p", col=NA)
    res.1 <- res[res$padj>=0.05 & abs(res$logFC)<=1, ]
    if(nrow(res.1)>0) points(x=res.1$logFC, y=res.1$qvalue, pch = 20, cex=.5, col="gray80")
    if(nrow(res.s)>0) points(x=res.s$logFC, y=res.s$qvalue, pch = 19, cex=.5, col=ifelse(res.s$logFC>0, "brown", "darkblue"))
    dev.off()
    res$qvalue <- -10*log10(res$PValue)
    png(fname(name, "png", "Volcano-plot", name))
    plot(x=res$logFC, y=res$qvalue,
         main = paste("Volcano plot for", name),
         xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
         type = "p", col=NA)
    res.1 <- res[res$padj>=0.05 & abs(res$logFC)<=1, ]
    if(nrow(res.1)>0) points(x=res.1$logFC, y=res.1$qvalue, pch = 20, cex=.5, col="gray80")
    if(nrow(res.s)>0) points(x=res.s$logFC, y=res.s$qvalue, pch = 19, cex=.5, col=ifelse(res.s$logFC>0, "brown", "darkblue"))
    dev.off()
  })
}




