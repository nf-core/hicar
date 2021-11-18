// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_CHIPPEAKANNO {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) }

    conda (params.enable_conda ? "bioconda::bioconductor-chippeakanno=3.26.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0"
    }

    input:
    tuple val(bin_size), path(diff)
    path gtf

    output:
    tuple val(bin_size), path("${prefix}/anno/*"), emit: anno
    tuple val(bin_size), path("${prefix}/anno/**.anno.csv"), emit: csv
    path "${prefix}/anno/*.png", optional:true, emit: png
    path "versions.yml"                       , emit: versions

    script:
    prefix   = options.suffix ? "${options.suffix}${bin_size}" : "diffhic_bin${bin_size}"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 call ChIPpeakAnno
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("ChIPpeakAnno", "rtracklayer", "GenomicFeatures", "ggplot2")
    versions <- c("${getProcessName(task.process)}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    gtf <- "${gtf}"
    pf <- file.path("${prefix}", "anno")
    bin_size <- "${prefix}"

    detbl <- dir(".", "DEtable.*.csv|sig3Dinteractions.bedpe|peaks",
                recursive = TRUE, full.names = TRUE)
    detbl <- detbl[!grepl("anno.csv", detbl)] ## in case of re-run

    txdb <- makeTxDbFromGFF(gtf) ## create annotation data from gtf file
    gtf <- import(gtf)
    id2symbol <- function(gtf){ ## convert entriz id to gene symbol
        if(is.null(gtf\$gene_name)) return(NULL)
        x <- data.frame(id=gtf\$gene_id, symbol=gtf\$gene_name)
        x <- unique(x)
        x <- x[!duplicated(x\$id), ]
        x <- x[!is.na(x\$id), , drop=FALSE]
        if(nrow(x)==0) return(NULL)
        y <- x\$symbol
        names(y) <- x\$id
        y
    }
    id2symbol <- id2symbol(gtf)
    anno <- toGRanges(txdb)
    promoters <- promoters(anno, upstream=2000, downstream=500)
    resList <- list() # save annotation results to a list
    peaks <- list()
    promoterList <- list() # list to save the distal sites of promoters interactions

    dir.create(pf, showWarnings = FALSE, recursive = TRUE)
    for(det in detbl){
        if(grepl("csv\$", det)) {
            DB <- read.csv(det)
        }else{
            if(grepl("peaks\$", det)){
                DB <- read.table(det, header=TRUE)
            }else{
                DB <- read.delim(det)
            }
        }
        if(nrow(DB)<1) next
        rownames(DB) <- paste0("p", seq.int(nrow(DB)))
        DB.gr1 <- with(DB, GRanges(chr1, IRanges(start1, end1, name=rownames(DB))))
        DB.gr2 <- with(DB, GRanges(chr2, IRanges(start2, end2, name=rownames(DB))))
        # Annotation
        DB.anno1 <- annotatePeakInBatch(DB.gr1, AnnotationData = anno,
                                        output = "both",
                                        PeakLocForDistance = "middle",
                                        FeatureLocForDistance = "TSS",
                                        ignore.strand = TRUE)
        if(length(id2symbol)>0) DB.anno1\$symbol[!is.na(DB.anno1\$feature)] <- id2symbol[DB.anno1\$feature[!is.na(DB.anno1\$feature)]]
        DB.anno2 <- annotatePeakInBatch(DB.gr2, AnnotationData = anno,
                                        output = "both",
                                        PeakLocForDistance = "middle",
                                        FeatureLocForDistance = "TSS",
                                        ignore.strand = TRUE)
        if(length(id2symbol)>0) DB.anno2\$symbol[!is.na(DB.anno2\$feature)] <- id2symbol[DB.anno2\$feature[!is.na(DB.anno2\$feature)]]
        groupName <- sub(".sig3Dinteractions.bedpe|csv", "", basename(det))
        if(grepl("padj", det)){
            resList[[groupName]] <- c(DB.anno1, DB.anno2)
        }else{
            peaks[[groupName]] <- unique(c(DB.gr1, DB.gr2))
            ol1 <- findOverlaps(DB.gr1, promoters)
            ol2 <- findOverlaps(DB.gr2, promoters)
            promoterList[[groupName]] <- unique(c(DB.gr2[unique(queryHits(ol1))], DB.gr1[unique(queryHits(ol2))]))
        }
        # Summary the annotations
        DB.anno1 <- mcols(DB.anno1)
        DB.anno2 <- mcols(DB.anno2)
        DB.anno <- merge(DB.anno1, DB.anno2, by="peak",
                        suffixes = c(".anchor1",".anchor2"))
        DB <- cbind(DB[DB.anno\$peak, ], DB.anno)
        pff <- file.path(pf, sub(".(csv|bedpe|peaks)", ".anno.csv", det))
        dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
        write.csv(DB, pff, row.names = FALSE)
    }


    if(packageVersion("ChIPpeakAnno")>="3.23.12"){
        if(length(resList)>0){
            if(is.list(resList)){
                resList <- GRangesList(resList[lengths(resList)>0])
            }
            out <- genomicElementDistribution(resList,
                                            TxDb = txdb,
                                            promoterRegion=c(upstream=2000, downstream=500),
                                            geneDownstream=c(upstream=0, downstream=2000),
                                            promoterLevel=list(
                                            # from 5' -> 3', fixed precedence 3' -> 5'
                                                breaks = c(-2000, -1000, -500, 0, 500),
                                                labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                        "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC", "#FFCA99",
                                                        "#FFAD65", "#FF8E32")),
                                            plot = FALSE)

            ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".pdf")), plot=out\$plot, width=9, height=9)
            ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".png")), plot=out\$plot)
            out <- metagenePlot(resList, txdb)
            ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".pdf")), plot=out, width=9, height=9)
            ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".png")), plot=out)
        }
        if(length(peaks)>0){
            peaks <- GRangesList(peaks[lengths(peaks)>0])
            out <- genomicElementDistribution(peaks,
                                            TxDb = txdb,
                                            promoterRegion=c(upstream=2000, downstream=500),
                                            geneDownstream=c(upstream=0, downstream=2000),
                                            promoterLevel=list(
                                                # from 5' -> 3', fixed precedence 3' -> 5'
                                                breaks = c(-2000, -1000, -500, 0, 500),
                                                labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                        "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC", "#FFCA99",
                                                        "#FFAD65", "#FF8E32")),
                                            plot = FALSE)

            ggsave(file.path(pf, paste0("genomicElementDistribuitonOfEachPeakList.", bin_size, ".pdf")), plot=out\$plot, width=9, height=9)
            ggsave(file.path(pf, paste0("genomicElementDistribuitonOfEachPeakList.", bin_size, ".png")), plot=out\$plot)

            out <- metagenePlot(peaks, txdb)
            ggsave(file.path(pf, paste0("metagenePlotToTSSOfEachPeakList.", bin_size, ".pdf")), plot=out, width=9, height=9)
            ggsave(file.path(pf, paste0("metagenePlotToTSSOfEachPeakList.", bin_size, ".png")), plot=out)

            if(length(peaks)<=5 && length(peaks)>1){
                ol <- findOverlapsOfPeaks(peaks)
                png(file.path(pf, paste0("vennDiagram.all.", bin_size, ".png")))
                vd <- makeVennDiagram(ol, connectedPeaks="keepAll")
                dev.off()
                write.csv(vd\$vennCounts, file.path(pf, paste0("vennDiagram.all.", bin_size, ".csv")), row.names=FALSE)
            }
        }
        if(length(promoterList)>0){
            promoterList <- GRangesList(promoterList[lengths(promoterList)>0])
            out <- genomicElementDistribution(promoterList,
                                            TxDb = txdb,
                                            promoterRegion=c(upstream=2000, downstream=500),
                                            geneDownstream=c(upstream=0, downstream=2000),
                                            promoterLevel=list(
                                                # from 5' -> 3', fixed precedence 3' -> 5'
                                                breaks = c(-2000, -1000, -500, 0, 500),
                                                labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                        "upstream <500b", "TSS - 500b"),
                                                colors = c("#FFE5CC", "#FFCA99",
                                                        "#FFAD65", "#FF8E32")),
                                            plot = FALSE)

            ggsave(file.path(pf, paste0("genomicElementDistribuitonOfremoteInteractionPeaks.", bin_size, ".pdf")), plot=out\$plot, width=9, height=9)
            ggsave(file.path(pf, paste0("genomicElementDistribuitonOfremoteInteractionPeaks.", bin_size, ".png")), plot=out\$plot)

            if(length(promoterList)<=5 && length(promoterList)>1){
                ol <- findOverlapsOfPeaks(promoterList)
                png(file.path(pf, paste0("vennDiagram.remote.interaction.peak.with.promoters.all.", bin_size, ".png")))
                vd <- makeVennDiagram(ol, connectedPeaks="keepAll")
                dev.off()
                write.csv(vd\$vennCounts, file.path(pf, paste0("vennDiagram.remote.interaction.peak.with.promoters.all.", bin_size, ".csv")), row.names=FALSE)
            }
        }
    }
    """
}
