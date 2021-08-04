#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on April. 29, 2021 call DESeq2
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(DESeq2)
writeLines(as.character(packageVersion("DESeq2")), "DESeq2.version.txt")

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
write.csv(sizeFactor, fname(NA, "csv", "size.factors"), row.names = TRUE)

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
  ## create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = cnts,
                                colData = coldata,
                                design = ~ condition)
  sizeFactors(dds) <- sizeFactor
  ## do differential analysis
  names(contrasts) <- vapply(contrasts, 
                             FUN=paste,
                             FUN.VALUE = character(1),
                             collapse = "-")
  # args$object <- dds
  # dds <- tryCatch(do.call(DESeq, args = args),
  #                 error = function(.e){
  #                   dds <- estimateSizeFactors(dds)
  #                   dds <- estimateDispersionsGeneEst(dds)
  #                   dispersions(dds) <- mcols(dds)$dispGeneEst
  #                   dds <- nbinomWaldTest(dds)
  #                 })
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  ## PCA
  pdf(fname(NA, "pdf", "PCA-plot"))
  plotPCA(rlog(dds, blind=FALSE),
          intgroup="condition", ntop = min(nrow(dds)-1, 500))
  dev.off()
  ## plot dispersion
  pdf(fname(NA, "pdf", "DispersionEstimate-plot"))
  plotDispEsts(dds)
  dev.off()
  res <- mapply(contrasts, names(contrasts), FUN = function(cont, name){
    rs <- results(dds, contrast = c("condition", cont))
    ## MA-plot
    pdf(fname(name, "pdf", "MA-plot", name))
    plotMA(rs)
    dev.off()
    ## pvale distribution
    pdf(fname(name, "pdf", "MA-plot", name))
    hist(rs$pvalue, breaks = 20)
    dev.off()
    ## save res
    res <- as.data.frame(rs)
    res <- cbind(peaks, res[names(peaks_id), ])
    write.csv(res, fname(name, "csv", "DESeq2.DEtable", name), row.names = FALSE)
    ## save metadata
    write.csv(elementMetadata(rs), fname(name, "csv", "DESeq2.metadata", name), row.names = FALSE)
    ## save subset results
    res.s <- res[res$padj<0.05 & abs(res$log2FoldChange)>1, ]
    write.csv(res.s, fname(name, "csv", "DESeq2.DEtable", name, "padj0.05.lfc1"), row.names = FALSE)
    ## Volcano plot
    res$qvalue <- -10*log10(res$pvalue)
    pdf(fname(name, "pdf", "Volcano-plot", name))
    plot(x=res$log2FoldChange, y=res$qvalue,
         main = paste("Volcano plot for", name),
         xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
         type = "p", col=NA)
    res.1 <- res[res$padj>=0.05 & abs(res$log2FoldChange)<=1, ]
    if(nrow(res.1)>0) points(x=res.1$log2FoldChange, y=res.1$qvalue, pch = 20, cex=.5, col="gray80")
    if(nrow(res.s)>0) points(x=res.s$log2FoldChange, y=res.s$qvalue, pch = 19, cex=.5, col=ifelse(res.s$log2FoldChange>0, "brown", "darkblue"))
    dev.off()
    res$qvalue <- -10*log10(res$pvalue)
    png(fname(name, "png", "Volcano-plot", name))
    plot(x=res$log2FoldChange, y=res$qvalue,
         main = paste("Volcano plot for", name),
         xlab = "log2 Fold Change", ylab = "-10*log10(P-value)",
         type = "p", col=NA)
    res.1 <- res[res$padj>=0.05 & abs(res$log2FoldChange)<=1, ]
    if(nrow(res.1)>0) points(x=res.1$log2FoldChange, y=res.1$qvalue, pch = 20, cex=.5, col="gray80")
    if(nrow(res.s)>0) points(x=res.s$log2FoldChange, y=res.s$qvalue, pch = 19, cex=.5, col=ifelse(res.s$log2FoldChange>0, "brown", "darkblue"))
    dev.off()
  })
}




