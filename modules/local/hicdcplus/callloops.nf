process HICDCPLUS_CALL_LOOPS {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'
    label 'error_ignore'

    conda "bioconda::bioconductor-hicdcplus=1.2.1"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-hicdcplus:1.2.1--r41h619a076_0' :
        'biocontainers/bioconductor-hicdcplus:1.2.1--r41h619a076_0' }"

    input:
    tuple val(meta), path(peaks), path(features), path(pairs), path(chromsize)

    output:
    tuple val(meta), val(meta.bin), path("*.interactions.bedpe")     , emit: interactions
    tuple val(meta), val(meta.bin), path("*_hicdcplus_result.txt.gz"), emit: full_results
    tuple val(meta), val(meta.bin), path("*.interactions.txt")       , emit: filtered_res
    path "versions.yml"                                              , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on July 13, 2022 for HiC-DC+ to call signficant interactions
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################

    pkgs <- c("HiCDCPlus", "dplyr", "InteractionSet", "GenomicRanges")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    Dthreshold <- 2e+06
    binsize <- ${meta.bin}
    NCORE <- ${task.cpus}
    allvalidpairs_paths <- strsplit("$pairs", "\\\\s+")[[1]]
    FDR <- 0.05
    HiCDCPlus_options <- NULL
    ## handle args
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
    option_list <- list("dthreshold"=c("--dthreshold", "-d", "numeric"),
                        "fdr"=c("--fdr", "-f", "numeric"),
                        "hicdcplus_options"=c("--hicdcplus_options", "-w", "character"))
    opt <- parse_args(option_list, args)
    if(!is.null(opt\$dthreshold)){
        Dthreshold <- opt\$dthreshold
    }
    if(!is.null(opt\$fdr)){
        FDR <- opt\$fdr
    }
    if(!is.null(opt\$hicdcplus_options)){
        HiCDCPlus_options <- opt\$hicdcplus_options
    }

    ## step 1, generate features
    bintolen <- data.table::fread("$features")
    bintolen <- bintolen %>%
        tidyr::separate(.data\$bins, c("chr", "start", "end"), "-") %>%
            dplyr::mutate(start = floor(as.numeric(.data\$start)/1000) *
                1000, end = as.numeric(.data\$end))
    seql <- read.delim("$chromsize", header=FALSE, row.names=1)
    seql <- GRanges(rownames(seql), IRanges(1, seql[, 1]))
    seql <- seql[seqnames(seql) %in% standardChromosomes(seql)]
    seql <- keepStandardChromosomes(seql)
    seqlen <- end(seql)
    names(seqlen) <- seqnames(seql)
    generate_binned_gi_list <- function(binsize, chrs = NULL, Dthreshold = 2e+06, seqlens) {
        gi_list <- list()
        for (chrom in chrs) {
            seqlen <- seqlens[chrom]
            numbins <- ceiling(seqlen/binsize)
            maxbins <- min(round(Dthreshold/binsize), numbins)
            all.regions <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = seq(0, (numbins - 1), 1) * binsize, end = pmin(seq(1, numbins, 1) * binsize, seqlen)))
            index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
            index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
            gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
            mcols(gi_list[[chrom]])\$D <- InteractionSet::pairdist(gi_list[[chrom]])
            gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])\$D <= Dthreshold]
        }
        return(gi_list)
    }
    gi_list <- generate_binned_gi_list(binsize, names(seqlen), Dthreshold, seqlen)
    ## filter the list to avoid the fitting error
    gi_list_D <- vapply(gi_list, function(.ele){
        sum(mcols(.ele)[["D"]]>0)>300 ## determined by ssize = 0.01, keep at least 3 point for model.
    }, logical(1L))
    gi_list <- gi_list[gi_list_D]
    # add counts
    ## re-arrange the pairs to a fake hic-pro validatedpair file
    for(allvalidpairs_path in allvalidpairs_paths){
        line1 <- data.table::fread(allvalidpairs_path, sep = "\\t", header = FALSE, stringsAsFactors = FALSE, nrows=1)
        if(line1[[1]]==line1[[4]]){
            allvalidpairs <- data.table::fread(allvalidpairs_path, sep = "\\t", header = FALSE, stringsAsFactors = FALSE)
            allvalidpairs <- allvalidpairs[rep(seq.int(nrow(allvalidpairs)), as.numeric(allvalidpairs\$V7)),  drop=FALSE]
            allvalidpairs <- cbind('.', allvalidpairs[, 1], rowMeans(allvalidpairs[, 2:3]), '.', allvalidpairs[, 4], rowMeans(allvalidpairs[, 5:6]))
            data.table::fwrite(allvalidpairs, 'tmp.txt', sep = '\\t', col.names = FALSE, row.names = FALSE, append = TRUE)
        }
    }
    gi_list <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path="tmp.txt")
    unlink("tmp.txt")
    #expand features for modeling
    gi_list <- expand_1D_features(gi_list)
    #run HiC-DC+
    set.seed(1) #HiC-DC downsamples rows for modeling
    if(!is.null(HiCDCPlus_options)){
        #HiCDCPlus_parallel runs in parallel across ncores
        gi_list<-HiCDCPlus_parallel(gi_list, ncore=NCORE, HiCDCPlus_options)
    }else{
        gi_list<-HiCDCPlus_parallel(gi_list, ncore=NCORE)
    }
    #write results to a text file for differential analysis
    gi_list_write(gi_list,fname="${prefix}_hicdcplus_result.txt.gz")
    #export the loops files
    gi_list <- Reduce(c, gi_list)
    gi_list <- gi_list[gi_list\$qvalue<FDR]
    write.table(as.data.frame(gi_list), "${prefix}.interactions.txt", sep="\\t", row.names=FALSE)
    mcols(gi_list) <- DataFrame(score = gi_list\$mu)
    rtracklayer::export(gi_list, "${prefix}.interactions.bedpe")
    """
}
