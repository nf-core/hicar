process ATACQC {
    label 'process_medium'

    conda "bioconda::bioconductor-atacseqqc=1.16.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0' :
        'biocontainers/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0' }"

    input:
    path peaks
    path beds
    path gtf

    output:
    path "*.csv"                   , emit: stats
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Sept. 10, 2021 stats for ATAC reads
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("rtracklayer", "GenomicFeatures", "ATACseqQC")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml
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
    option_list <- list("pattern"=c("--pattern", "-p", "character"))
    opt <- parse_args(option_list, args)
    pattern <- opt[['pattern']] #"merged.ATAC.bed.gz"
    gtf <- "$gtf"

    readsFiles <- dir(".", pattern) ## postfix from mergedreads.nf
    peaksFiles <- dir(".", "narrowPeak|broadPeak") ## output from macs2

    names(readsFiles) <- sub(pattern, "", readsFiles)
    names(peaksFiles) <- sub("_peaks.*Peak", "", peaksFiles)

    N <- intersect(names(readsFiles), names(peaksFiles))
    if(length(N)==0){
        if(length(peaksFiles)==1){ # for user defined peaks.
            peaksFiles <- rep(peaksFiles, length(readsFiles))
            names(peaksFiles) <- names(readsFiles)
            N <- intersect(names(readsFiles), names(peaksFiles))
        }
        stopifnot("no peak and signal pairs"=length(N)>0)
    }

    peaksFiles <- peaksFiles[N]
    readsFiles <- readsFiles[N]

    ## import reads
    readls <- lapply(readsFiles, function(f){
        reads <- read.table(f, colClasses=c("character", "integer", "NULL",
                                "NULL", "NULL", "character"))
        reads <- GRanges(reads[, 1],
                        IRanges(as.numeric(reads[, 2])+1, as.numeric(reads[, 2])+150),
                        strand = reads[, 3])
    })

    ## import peaks
    peakls <- lapply(peaksFiles, import)

    txdb <- makeTxDbFromGFF(gtf) #for TSS
    txs <- exons(txdb)

    stats <- mapply(peakls, readls, FUN = function(peaks, reads){
        ## calculate FRiP score (fraction of reads in peaks), must over 1%
        readsInPeaks <- countOverlaps(peaks, reads)
        FRiP <- 100*sum(readsInPeaks)/length(reads)

        reads <- as(reads, "GAlignments")
        ## calculate Transcription Start Site Enrichment Score
        tsse <- TSSEscore(reads, txs)

        ## Promoter/Transcript body (PT) score
        pt <- PTscore(reads, txs)
        pt <- pt[!is.na(pt\$promoter) & !is.na(pt\$transcriptBody) & !is.na(pt\$PT_score)]
        pt <- pt[(pt\$promoter>0 | pt\$transcriptBody>0) & pt\$PT_score!=0]
        promoterEnriched <- table(pt\$PT_score>0)
        names(promoterEnriched) <-
            c("FALSE"="bodyEnrich", "TRUE"="promoterEnrich")[names(promoterEnriched)]
        promoterEnriched <-
            c(promoterEnriched,
                prop.test=prop.test(cbind(table(pt\$PT_score>0), c(50, 50)))\$p.value)
        list(tsse_FRiP = c(TSSEscore=tsse\$TSSEscore, FRiP=FRiP, promoterEnriched),
            tsseValues = tsse\$values)
    }, SIMPLIFY = FALSE)

    ## FRiP, TSSEscore table
    tsse_FRiP <- do.call(rbind, lapply(stats, function(.ele) .ele\$tsse_FRiP))
    write.csv(tsse_FRiP, "TSSEscore_FRiP.csv")

    ## for TSSEscore to TSS plots
    tsse <- do.call(rbind, lapply(stats, function(.ele) .ele\$tsseValues))
    if(ncol(tsse)==20){
        colnames(tsse) <- 100*(-9:10-.5)
        write.csv(tsse, "aggregateTSSEscoreToTSS.csv")
    }
    """
}
