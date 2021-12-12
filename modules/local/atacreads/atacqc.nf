// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process ATACQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-atacseqqc=1.16.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0"
    }

    input:
    path peaks
    path beds
    path gtf

    output:
    path "*.csv"                   , emit: stats
    path "versions.yml"            , emit: versions

    script:
    def software = "ATACseqQC"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Sept. 10, 2021 stats for ATAC reads
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("rtracklayer", "GenomicFeatures", "ATACseqQC")
    versions <- c("${getProcessName(task.process)}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    gtf <- "$gtf"

    readsFiles <- dir(".", "merged.ATAC.bed.gz") ## postfix from mergedreads.nf
    peaksFiles <- dir(".", "narrowPeak|broadPeak") ## output from macs2

    names(readsFiles) <- sub(".merged.ATAC.bed.gz", "", readsFiles)
    names(peaksFiles) <- sub("_peaks.*Peak", "", peaksFiles)

    N <- intersect(names(readsFiles), names(peaksFiles))
    stopifnot("no peak and signal pairs"=length(N)>0)

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
