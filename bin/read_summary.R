#!/usr/bin/env Rscript

fs <- dir(".", "*.reads_stats.csv")
df <- do.call(rbind, lapply(fs, read.csv))
con <- file("reads_summary.csv", open="wt")
write.csv(df, con, row.names=FALSE)
close(con)
