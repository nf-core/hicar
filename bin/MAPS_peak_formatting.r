#!/usr/bin/env Rscript

#########################################
# Author: [Ivan Juric](https://github.com/ijuric)
# File: MAPS_peak_formatting.r
# Source: https://github.com/ijuric/MAPS/blob/master/bin/MAPS/MAPS_peak_formatting.r
# Data: 11/08/2021
# modified by Jianhong:
# ## 1. set the Rscript environment.
# ## 2. simplify the input and output parameters by input filename and output filename
# ## 3. handle multiple input files
# ## 4. handle the empty input file
#########################################

options("scipen"=999)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

RESOLUTION = as.numeric(args[1])
infs = args[-1]

for(inf in infs){
    peaks_raw = read.table(inf,header=T, stringsAsFactors=F)

    peaks = as.data.table(subset(peaks_raw, ClusterType != 'Singleton' | (ClusterType == 'Singleton' & fdr < 1e-4))) # remove singletons
    if(nrow(peaks)==0){
        peaks$bin1_end <- peaks$bin2_end <- peaks$summit <- numeric(0)
    }else{
        peaks[, summit := 1*(fdr == min(fdr)), by = lab]
        peaks$summit[ peaks$ClusterType == 'Singleton'] = 1

        singleton_labs = peaks$lab[ peaks$ClusterType == 'Singleton']
        peaks$lab[ peaks$ClusterType == 'Singleton'] = paste(singleton_labs, 1:length(singleton_labs),sep='')

        peaks$bin1_end = peaks$bin1_mid + RESOLUTION
        peaks$bin2_end = peaks$bin2_mid + RESOLUTION
    }
    peaks_final = subset(peaks, select = c("chr", "bin1_mid", "bin1_end", "chr", "bin2_mid", "bin2_end", "count", "expected2", "fdr", "lab", "ClusterSize", "ClusterType", "NegLog10P", "summit"))
    colnames(peaks_final) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'count', 'expected', 'fdr', 'ClusterLabel', 'ClusterSize', 'ClusterType', 'ClusterNegLog10P', 'ClusterSummit')
    fout = sub(".peaks",'.sig3Dinteractions.bedpe',inf)
    write.table(peaks_final, fout, row.names = FALSE, col.names = TRUE, quote=FALSE, sep='\t')
}
