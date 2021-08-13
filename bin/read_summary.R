#!/usr/bin/env Rscript

fs <- dir(".", "*.reads_stats.csv")
df <- do.call(rbind, lapply(fs, read.csv))
con <- file("reads_summary.csv", open="wt")
write.csv(df, con, row.names=FALSE)
close(con)

fs <- dir(".", "*.summary.out$")
df <- t(do.call(cbind, lapply(fs, read.delim, header=FALSE, row.names=1)))
df <- sub(",", "", df)
mode(df) <- "numeric"
df <- cbind(sample=sub(".summary.out", "", basename(fs)), df)
con <- file("pairsqc_summary_out.csv", open="wt")
write.csv(df, con, row.names=FALSE)
close(con)
