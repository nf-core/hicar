#!/usr/bin/env Rscript
## generate the statistis for each samples
args <- commandArgs(trailingOnly=TRUE)
raw = args[1]
dedup = args[2]
out = args[3]

sample_name <- sub(".raw.pairsam.stat", "", basename(raw))
getDat <- function(f){
    dat <- read.delim(f, header=FALSE)
    res <- dat[, 2]
    names(res) <- dat[, 1]
    return(res)
}
all_pairs <- getDat(raw)
dep_pairs <- getDat(dedup)

df <- data.frame(sample=sample_name,
        total=all_pairs["total"],
        duplicate=dep_pairs['total_dups'],
        non_duplicated=dep_pairs['total_nodups'],
        trans=dep_pairs['trans'],
        cis=dep_pairs['cis'],
        longRange=dep_pairs['cis_20kb+'])
write.csv(df, out, row.names=FALSE)
