process ASSIGN_TYPE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bioconductor-chippeakanno=3.26.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0' :
        'biocontainers/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("summary.*.txt"), optional: true, emit: summary
    tuple val(meta), path("*.peaks")      , optional: true, emit: peak
    tuple val(meta), path("*.bedpe")      , optional: true, emit: bedpe
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 assign interacion type for the peaks
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("graph", "RBGL", "InteractionSet")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    ## Options
    ## make_option(c("-c", "--count_cutoff"), type="integer", default=12, help="count cutoff, default 12", metavar="integer")
    ## make_option(c("-r", "--ratio_cutoff"), type="numeric", default=2.0, help="ratio cutoff, default 2.0", metavar="float")
    ## make_option(c("-f", "--fdr"), type="integer", default=2, help="-log10(fdr) cutoff, default 2", metavar="integer")
    ## make_option(c("-i", "--interactions"), type="character", default=NULL, help="interactions output by call hipeak", metavar="string")
    ## make_option(c("-o", "--output"), type="character", default="peaks", help="sample name of the output prefix", metavar="string")

    OUTPUT = "."
    GROUP_ID = "$meta.id"
    COUNT_CUTOFF = 3
    RATIO_CUTOFF = 2.0
    FDR = 2
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
    option_list <- list("count_cutoff"=c("--count_cutoff", "-c", "integer"),
                        "ratio_cutoff"=c("--ratio_cutoff", "-r", "numeric"),
                        "fdr"=c("--fdr", "-f", "integer"))
    opt <- parse_args(option_list, strsplit("$args", "\\\\s+")[[1]])
    if(!is.null(opt\$count_cutoff)){
        COUNT_CUTOFF <- opt\$count_cutoff
    }
    if(!is.null(opt\$ratio_cutoff)){
        RATIO_CUTOFF <- opt\$ratio_cutoff
    }
    if(!is.null(opt\$fdr)){
        FDR <- opt\$fdr
    }

    peaks <- read.csv("$counts")

    if(!all(c("chr1", "start1", "end1", "width1",
            "chr2", "start2", "end2", 'width2',
            "count", "logl", "logn", "loggc", "logm", "logdist", 'logShortCount',
            "ratio2", 'fdr') %in% colnames(peaks))){
        stop("count table is not in correct format.")
    }

    classify_peaks <- function(final) {
        # group the interactions
        gi <- with(final, GInteractions(GRanges(chr1, IRanges(start1, end1)), GRanges(chr2, IRanges(start2, end2))))
        ol1 <- findOverlaps(first(gi), drop.self = TRUE, drop.redundant = TRUE)
        ol2 <- findOverlaps(second(gi), drop.self = TRUE, drop.redundant = TRUE)
        ol <- unique(c(queryHits(ol1), subjectHits(ol1), queryHits(ol2), subjectHits(ol2)))
        ol_ <- seq_along(gi)[-ol]

        group <- unique(rbind(as.data.frame(ol1), as.data.frame(ol2)))
        colnames(group) <- c("from", "to")
        group\$weight <- rep(1L, nrow(group))
        group <- split(group, seqnames(first(gi)[group\$from]))
        group <- mapply(group, names(group), FUN=function(.ele, .name){
            if(length(unique(c(.ele[, 1], .ele[, 2])))>sqrt(.Machine\$integer.max)){
                .nodes <- unique(c(.ele[, 1], .ele[, 2]))
                .max_nodes <- floor(sqrt(.Machine\$integer.max)/2)
                .nodes_group <- rep(seq.int(ceiling(length(.nodes)/.max_nodes)),
                                    each=.max_nodes)[seq_along(.nodes)]
                names(.nodes_group) <- .nodes
                .ele <- split(.ele, paste(.nodes_group[as.character(.ele[, 1])],
                                        .nodes_group[as.character(.ele[, 2])],
                                        sep="_"))
                .ele <- lapply(.ele, function(.e){
                    .e <- graphBAM(.e)
                    .e <- connectedComp(ugraph(.e))
                    .e <- lapply(.e, as.numeric)
                })
                .g <- mapply(.ele, seq_along(.ele), FUN=function(.e, .id){
                    data.frame( nodes=unlist(.e),
                                group=rep(paste(.id, names(.e), sep="_"), lengths(.e)))
                }, SIMPLIFY=FALSE)
                .g <- do.call(rbind, .g)
                .gs <- split(.g[, 2], .g[, 1])
                .gs <- lapply(.gs, unique)
                .gs <- .gs[!duplicated(.gs)]
                while(any(lengths(.gs)>1)){
                    ## merge parents
                    .gsn <- vapply(.gs, FUN=function(.e) sort(.e)[1], FUN.VALUE=character(1))
                    .gsn <- rep(.gsn, lengths(.gs))
                    names(.gsn) <- unlist(.gs)
                    .gsn <- .gsn[names(.gsn)!=.gsn]
                    .k <- .g[, "group"] %in% names(.gsn)
                    .g[.k, "group"] <- .gsn[.g[.k, "group"]]
                    .gs <- split(.g[, 2], .g[, 1])
                    .gs <- lapply(.gs, unique)
                    .gs <- .gs[!duplicated(.gs)]
                }
                .g <- unique(.g)
                .ele <- split(.g[, "nodes"], .g[, "group"])
                names(.ele) <- seq_along(.ele)
                rm(.g, .gs, .gsn, .nodes_group, .nodes)
            }else{
                .ele <- graphBAM(.ele)
                .ele <- connectedComp(ugraph(.ele))
                .ele <- lapply(.ele, as.numeric)
            }
            data.frame( id=unlist(.ele),
                        g=rep(paste(.name, seq_along(.ele), sep="_"), lengths(.ele)))
        }, SIMPLIFY=FALSE)
        group <- do.call(rbind, group)

        final\$Cluster <- NA
        final\$Cluster[group\$id] <- group\$g
        final\$ClusterSize <- 0
        final\$ClusterSize[group\$id] <- table(group\$g)[group\$g]
        if(any(is.na(final\$Cluster))) final\$Cluster[is.na(final\$Cluster)] <- paste0("Singleton_", seq.int(sum(is.na(final\$Cluster))))
        final\$NegLog10P <- -log10( final\$p_val_reg2 )
        final\$NegLog10P[is.na(final\$NegLog10P)] <- 0
        final\$NegLog10P[is.infinite(final\$NegLog10P)] <- max(final\$NegLog10P[!is.infinite(final\$NegLog10P)]+1)
        NegLog10P <- rowsum(final\$NegLog10P, final\$Cluster)
        final\$NegLog10P <- NegLog10P[final\$Cluster, 1]

        x <- unique( final[ final\$ClusterSize != 0, c('chr1', 'Cluster', 'NegLog10P', 'ClusterSize')] )
        if(nrow(x)==0){
            final\$ClusterType <- 'Singleton'
            return(final)
        }

        # sort rows by cumulative -log10 P-value
        x <- x[ order(x\$NegLog10P) ,]
        y<-sort(x\$NegLog10P)
        z<-cbind( seq(1,length(y),1), y )

        # keep a record of z before normalization
        z0 <- z

        z[,1]<-z[,1]/max(z[,1], na.rm=TRUE)
        z[,2]<-z[,2]/max(z[,2], na.rm=TRUE)

        u<-z
        u[,1] <-  1/sqrt(2)*z[,1] + 1/sqrt(2)*z[,2]
        u[,2] <- -1/sqrt(2)*z[,1] + 1/sqrt(2)*z[,2]

        v<-cbind(u, seq(1,nrow(u),1) )
        RefPoint <- v[ v[,2]==min(v[,2], na.rm=TRUE) , 3]
        RefValue <- z0[RefPoint,2]

        # define peak cluster type
        final\$ClusterType <- rep(NA, nrow(final))
        if(length(ol_)) final\$ClusterType[ ol_ ] <- 'Singleton'
        if(length(ol)){
            final\$ClusterType[ seq_along(gi) %in% ol & final\$NegLog10P < RefValue ] <-  'SharpPeak'
            final\$ClusterType[ seq_along(gi) %in% ol & final\$NegLog10P >= RefValue ] <- 'BroadPeak'
        }
        return(final)
    }

    peaks <- if(nrow(peaks)>0) subset(peaks, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > FDR) else data.frame()
    if (dim(peaks)[1] == 0) {
        print(paste('ERROR caller_hipeak.r: 0 bin pairs with count >= ',COUNT_CUTOFF,' observed/expected ratio >= ',RATIO_CUTOFF,' and -log10(fdr) > ',FDR,sep=''))
        quit()
    }

    peaks = classify_peaks(peaks)

    outf_name = paste(GROUP_ID, '.',FDR,'.peaks',sep='')
    dir.create(OUTPUT, recursive=TRUE)
    peaks <- unique(peaks)
    write.table(peaks, file.path(OUTPUT, outf_name),
                row.names = FALSE, col.names = TRUE, quote=FALSE)
    peaks1 <- cbind(peaks[, c("chr1", "start1", "end1", "chr2", "start2", "end2")], "*", peaks[, "NegLog10P", drop=FALSE])
    peaks1 <- unique(peaks1)
    write.table(peaks1,
                file.path(OUTPUT, paste0(GROUP_ID, '.', FDR, '.bedpe')),
                row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")

    summary_all_runs <- split(peaks, peaks\$ClusterType)
    summary_all_runs <- lapply(summary_all_runs, function(.ele){
        c(count = nrow(.ele),
        minWidth1 = min(.ele\$width1),
        medianWidth1 = median(.ele\$width1),
        maxWidth1 = max(.ele\$width1),
        minWidth2 = min(.ele\$width2),
        medianWidth2 = median(.ele\$width2),
        maxWidth2 = max(.ele\$width2),
        minFoldChange = min(.ele\$ratio2),
        medianFoldChange = median(.ele\$ratio2),
        maxFoldChange = max(.ele\$ratio2))
    })
    summary_all_runs <- do.call(rbind, summary_all_runs)
    summary_outf_name = paste('summary.',GROUP_ID,'.txt',sep='')
    write.table(summary_all_runs, file.path(OUTPUT, summary_outf_name), row.names = TRUE, col.names = TRUE, quote=FALSE)
    """
}
