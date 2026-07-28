process BIOC_CHIPPEAKANNO {
    tag "$foldername"
    label 'process_medium'
    //label 'error_ignore'

    conda "bioconda::bioconductor-chippeakanno=3.32.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.32.0--r42hdfd78af_0' :
        'biocontainers/bioconductor-chippeakanno:3.32.0--r42hdfd78af_0' }"

    input:
    tuple val(foldername), path(diff), path(peak)
    path gtf
    val maps_3d_ext

    output:
    tuple val(foldername), path("${prefix}/*")          , emit: anno
    tuple val(foldername), path("${prefix}/**.anno.csv"), emit: csv
    path "${prefix}/*.png", optional:true               , emit: png
    path "versions.yml"                                 , emit: versions

    script:
    prefix   = task.ext.prefix ?: "$foldername"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 call ChIPpeakAnno
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("ChIPpeakAnno", "rtracklayer", "GenomicFeatures", "ggplot2")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    gtf <- "${gtf}"
    pf <- file.path("${prefix}")
    bin_size <- basename("${prefix}")
    anchor_peak <- "${peak}"

    detbl <- dir(".", "DEtable.*.csv|${maps_3d_ext}|peaks|bedpe|bed",
                recursive = TRUE, full.names = TRUE)
    detbl <- detbl[!grepl("anno.csv", detbl)] ## in case of re-run
    detbl <- detbl[!basename(detbl) %in% basename(anchor_peak)]
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
    R2_peaks <- list() # save annotation results of R2 peaks involved in the loops

    dir.create(pf, showWarnings = FALSE, recursive = TRUE)
    if(file.exists(anchor_peak)){
        R2 <- import(anchor_peak)
    }else{
        R2 <- GRanges()
    }
    for(det in detbl){
        DB <- data.frame()
        if(grepl("csv\$", det)) {
            DB <- read.csv(det)
        }else{
            if(grepl("peaks\$", det)){
                DB <- read.table(det, header=TRUE)
            }else{
                header <- tryCatch(
                    read.table(det, header=FALSE, nrow=1),
                    error = function(.e){
                        message(.e)
                    },
                    finally = NULL
                )
                if(length(header)>0){
                    hasHeader <- all(c("chr1", "start1", "end1",
                                        "chr2", "start2", "end2") %in%
                                        header[1, , drop=TRUE])
                    DB <- read.table(det, header = hasHeader,
                                        stringsAsFactors = FALSE)
                    if(!hasHeader){
                        colnames(DB)[1:6] <- c("chr1", "start1", "end1",
                                            "chr2", "start2", "end2")
                    }
                }
            }
        }
        if(nrow(DB)<1) next
        rownames(DB) <- paste0("p", seq.int(nrow(DB)))
        DB.gr1 <- with(DB, GRanges(chr1, IRanges(start1, end1, name=rownames(DB))))
        DB.gr2 <- with(DB, GRanges(chr2, IRanges(start2, end2, name=rownames(DB))))
        groupName <- gsub(".(sig3Dinteractions.pe.txt|csv|bedpe|txt)", "", basename(det))
        # Annotation
        if(length(R2)){
            ssrR2 <- seqlevelsStyle(R2)
            ssrDB <- seqlevelsStyle(DB.gr1)
            if(length(intersect(ssrR2, ssrDB))==0){
                seqlevelsStyle(R2) <- ssrDB[1]
            }
            ol1 <- findOverlaps(R2, DB.gr1)
            DB.gr1.R2 <- R2[queryHits(ol1)]
            DB.gr1.R2\$idx <- subjectHits(ol1)
            ol2 <- findOverlaps(R2, DB.gr2)
            DB.gr2.R2 <- R2[queryHits(ol2)]
            DB.gr2.R2\$idx <- subjectHits(ol2)
            R2_peaks[[groupName]] <- unique(c(DB.gr1.R2, DB.gr2.R2))
            DB.anno1 <- annotatePeakInBatch(DB.gr1.R2, AnnotationData = anno,
                                            output = "both",
                                            PeakLocForDistance = "middle",
                                            FeatureLocForDistance = "TSS",
                                            ignore.strand = TRUE)
            DB.anno2 <- annotatePeakInBatch(DB.gr2.R2, AnnotationData = anno,
                                            output = "both",
                                            PeakLocForDistance = "middle",
                                            FeatureLocForDistance = "TSS",
                                            ignore.strand = TRUE)
            m.DB.anno1 <- as.data.frame(DB.anno1)
            DB.anno1 <- DB.gr1[DB.anno1\$idx]
            DB.anno1\$peak <- names(DB.anno1)
            m.DB.anno1 <- m.DB.anno1[,!colnames(m.DB.anno1) %in% c("width", "strand", "name", "idx", "peak"), drop=FALSE]
            colnames(m.DB.anno1)[c(1, 2, 3)] <- paste0("ATAC_peak_", colnames(m.DB.anno1)[c(1, 2, 3)])
            mcols(DB.anno1) <- cbind(mcols(DB.anno1), m.DB.anno1)
            m.DB.anno2 <- as.data.frame(DB.anno2)
            DB.anno2 <- DB.gr2[DB.anno2\$idx]
            DB.anno2\$peak <- names(DB.anno2)
            m.DB.anno2 <- m.DB.anno2[,!colnames(m.DB.anno2) %in% c("width", "strand", "name", "idx", "peak"), drop=FALSE]
            colnames(m.DB.anno2)[c(1, 2, 3)] <- paste0("ATAC_peak_", colnames(m.DB.anno2)[c(1, 2, 3)])
            mcols(DB.anno2) <- cbind(mcols(DB.anno2), m.DB.anno2)
        }else{
            DB.anno1 <- annotatePeakInBatch(DB.gr1, AnnotationData = anno,
                                            output = "both",
                                            PeakLocForDistance = "middle",
                                            FeatureLocForDistance = "TSS",
                                            ignore.strand = TRUE)
            DB.anno2 <- annotatePeakInBatch(DB.gr2, AnnotationData = anno,
                                            output = "both",
                                            PeakLocForDistance = "middle",
                                            FeatureLocForDistance = "TSS",
                                            ignore.strand = TRUE)
        }
        ## unique annotation plots
        ## only take the nearest annotation
        DB.anno1.srt <- DB.anno1[order(abs(DB.anno1\$distancetoFeature))]
        DB.anno2.srt <- DB.anno2[order(abs(DB.anno2\$distancetoFeature))]
        DB.anno1.srt <- unique(DB.anno1.srt)
        DB.anno2.srt <- unique(DB.anno2.srt)
        DB.anno1.srt <- DB.anno1.srt[order(as.numeric(sub('p', '', names(DB.anno1.srt))))]
        DB.anno2.srt <- DB.anno2.srt[order(as.numeric(sub('p', '', names(DB.anno2.srt))))]
        intersect_names <- intersect(names(DB.anno1.srt), names(DB.anno2.srt))
        DB.anno1.srt <- DB.anno1.srt[intersect_names]
        DB.anno2.srt <- DB.anno2.srt[intersect_names]
        pff <- file.path(pf, sub(".(csv|bedpe|peaks|txt)\$", "", det))
        dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
        plotdata <- data.frame(
            distance1=DB.anno1.srt\$distancetoFeature,
            distance2=DB.anno2.srt\$distancetoFeature)
        plotdata\$r <- sqrt(plotdata\$distance1^2+plotdata\$distance2^2)
        plotdata\$type1 <- ifelse(abs(plotdata\$distance1) <= 2000, "promoter (2K)", "remote element")
        plotdata\$type2 <- ifelse(abs(plotdata\$distance2) <= 2000, "promoter (2K)", "remote element")
        type12 <- apply(plotdata[, c("type1", "type2")], 1, sort, simplify = FALSE)
        plotdata\$ChromosomeRegion <- vapply(type12, FUN=function(.ele) paste(.ele, collapse=" / "), FUN.VALUE = character(1L))
        log10transform <- function(x){
            ifelse(x==0, 0, sign(x)*log10(abs(x)))
        }
        plotdata\$log10_transformed_distance1 <- log10transform(plotdata\$distance1)
        plotdata\$log10_transformed_distance2 <- log10transform(plotdata\$distance2)
        plot <- ggplot(plotdata, aes(x=log10_transformed_distance1, y=log10_transformed_distance2, color=ChromosomeRegion)) +
            geom_point() + coord_fixed() + theme_bw()
        ggsave(paste0(pff, ".anno.distanceToTSS.xyscatter.pdf"), plot=plot, width=7, height=7)
        plot <- ggplot(plotdata, aes(x=ChromosomeRegion, fill=ChromosomeRegion)) +
            geom_bar(stat="count") + theme_bw() + theme(axis.text.x = element_blank())
        ggsave(paste0(pff, ".anno.distanceToTSS.barplot.pdf"), plot=plot, width=5, height=4)
        if(length(id2symbol)>0) DB.anno1.srt\$symbol[!is.na(DB.anno1.srt\$feature)] <- id2symbol[DB.anno1.srt\$feature[!is.na(DB.anno1.srt\$feature)]]
        if(length(id2symbol)>0) DB.anno2.srt\$symbol[!is.na(DB.anno2.srt\$feature)] <- id2symbol[DB.anno2.srt\$feature[!is.na(DB.anno2.srt\$feature)]]
        DB.anno1.srt <- mcols(DB.anno1.srt)
        DB.anno2.srt <- mcols(DB.anno2.srt)
        DB.anno.srt <- merge(DB.anno1.srt, DB.anno2.srt, by="peak",
            suffixes = c(".anchor1",".anchor2"),
            all = TRUE)
        DB.anno.srt <- cbind(DB[DB.anno.srt\$peak, , drop=FALSE], DB.anno.srt)
        DB.anno.srt <- DB.anno.srt[order(as.numeric(sub('p', '', DB.anno.srt\$peak))), , drop=FALSE]
        pff <- file.path(pf, sub(".(csv|bedpe|peaks|txt)\$", ".unique.anno.csv", det))
        write.csv(DB.anno.srt, pff, row.names = FALSE)

        if(length(id2symbol)>0) DB.anno1\$symbol[!is.na(DB.anno1\$feature)] <- id2symbol[DB.anno1\$feature[!is.na(DB.anno1\$feature)]]
        if(length(id2symbol)>0) DB.anno2\$symbol[!is.na(DB.anno2\$feature)] <- id2symbol[DB.anno2\$feature[!is.na(DB.anno2\$feature)]]
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
                        suffixes = c(".anchor1",".anchor2"),
                        all = TRUE)
        DB <- cbind(DB[DB.anno\$peak, ], DB.anno)
        pff <- file.path(pf, sub(".(csv|bedpe|peaks|txt)\$", ".anno.csv", det))
        dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
        write.csv(DB, pff, row.names = FALSE)
    }


    if(packageVersion("ChIPpeakAnno")>="3.23.12"){
        if(length(resList)>0){
            if(is.list(resList)){
                resList <- GRangesList(resList[lengths(resList)>0])
            }
            if(length(resList)>0){
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
                saveRDS(out\$peaks, file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".RDS")))
                ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".pdf")), plot=out\$plot, width=9, height=9)
                ggsave(file.path(pf, paste0("genomicElementDistribuiton.", bin_size, ".png")), plot=out\$plot)
                out <- metagenePlot(resList, txdb)
                ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".pdf")), plot=out, width=9, height=9)
                ggsave(file.path(pf, paste0("metagenePlotToTSS.", bin_size, ".png")), plot=out)
            }
        }
        if(length(peaks)>0){
            peaks <- GRangesList(peaks[lengths(peaks)>0])
            if(length(peaks)>0){
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
                saveRDS(out\$peaks, file.path(pf, paste0("genomicElementDistribuitonOfEachPeakList.", bin_size, ".RDS")))
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
        }
        if(length(promoterList)>0){
            promoterList <- GRangesList(promoterList[lengths(promoterList)>0])
            if(length(promoterList)>0){
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
                saveRDS(out\$peaks, file.path(pf, paste0("genomicElementDistribuitonOfremoteInteractionPeaks.", bin_size, ".RDS")))
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
        if(length(R2_peaks)>0){
            R2_peaks <- GRangesList(R2_peaks[lengths(R2_peaks)>0])
            if(length(R2_peaks)>0){
                out <- genomicElementDistribution(R2_peaks,
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
                saveRDS(out\$peaks, file.path(pf, paste0("genomicElementDistribuitonOfATACpeakInvolvedInInteractionPeaks.", bin_size, ".RDS")))
                ggsave(file.path(pf, paste0("genomicElementDistribuitonOfATACpeakInvolvedInInteractionPeaks.", bin_size, ".pdf")), plot=out\$plot, width=9, height=9)
                ggsave(file.path(pf, paste0("genomicElementDistribuitonOfATACpeakInvolvedInInteractionPeaks.", bin_size, ".png")), plot=out\$plot)

                if(length(R2_peaks)<=5 && length(R2_peaks)>1){
                    ol <- findOverlapsOfPeaks(R2_peaks)
                    png(file.path(pf, paste0("vennDiagram.ATACpeakInvolvedInInteractionPeaks.all.", bin_size, ".png")))
                    vd <- makeVennDiagram(ol, connectedPeaks="keepAll")
                    dev.off()
                    write.csv(vd\$vennCounts, file.path(pf, paste0("vennDiagram.ATACpeakInvolvedInInteractionPeaks.all.", bin_size, ".csv")), row.names=FALSE)
                }
            }
        }
    }
    """
}
