#!/usr/bin/env Rscript
## run example:
##  Rscript regression_and_peak_caller.r /home/jurici/work/PLACseq/MAPS_pipe/results/mESC_test/ MY_115.5k 5000 1 None pospoisson NA
##
## arguments:
## OUTPUT - output folder
## R2PEAK - R2 peak path
## R1PEAK - R1 peak path
## PAIRS - distal pairs path
## FASTA - genome fasta file path
## MAPPABILITY - mappability file path
## CUT - enzyme digest positions
## COUNT_CUTOFF - count cutoff, default 12
## RATIO_CUTOFF - ratio cutoff, default 2.0
## FDR - -log10(fdr) cutoff, default 2
## FILTER - file containing bins that need to be filtered out. Format: two columns "chrom", "bin". "chrom" contains 'chr1','chr2',.. "bin" is bin label
## regresison_type - pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson

library(VGAM)
library(MASS)
options(warn=-1)

### constants
COUNT_CUTOFF = 12
RATIO_CUTOFF = 2.0
FDR = 2

REG_TYPE = 'pospoisson'
###

args <- commandArgs(trailingOnly=TRUE)
fltr = data.frame(chr='chrNONE',bin=-1)

if (length(args) < 7 || length(args) > 12) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): OUTPUT, R2PEAK, R1PEAK, PAIRS, FASTA, MAPPABILITY, CUT, COUNT_CUTOFF, RATIO_CUTOFF, FDR, FILTER, regression_type.')
    print('FILTER is optional argument. Omitt it if no filtering required.')
    print(paste('Number of arguments entered:',length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    OUTPUT = args[1]
    R2PEAK = args[2]
    R1PEAK <- args[3]
    PAIRS <- args[4]
    FASTA <- args[5]
    MAPPABILITY <- args[6]
    CUT <- args[7]
    if(length(args)>7){
        COUNT_CUTOFF = as.numeric(args[8])
        RATIO_CUTOFF = as.numeric(args[9])
        FDR = as.numeric(args[10])
        print('filter used (if any):')
        if (length(args) > 10) {
            if (args[12] != 'None') {
                FILTER = args[11]
                fltr = read.table(FILTER,header=T)
                print(fltr)
            } else {
                print('None')
            }
            ## this done so that the script is compatible with previous run_pipeline scripts
            if (length(args) > 11) {
                if (args[12] != 'pospoisson' && args[12] != 'negbinom') {
                    print(paste('wrong regression choice. Your choice:', args[12], '. Avaiable choices: pospoisson or negbinom'),sep = ' ')
                    quit()
                }
                REG_TYPE = args[12]
            }
        } else {
            print('None')
        }
    }
}

## loading data
library(rtracklayer)
library(InteractionSet)
library(Biostrings)
library(Rsamtools)
library(graph)
library(RBGL)
pairs <- dir(PAIRS, full.names=TRUE)
R1PEAK <- import(R1PEAK)
R2PEAK <- import(R2PEAK)
## split by chromsome
R1PEAK <- split(R1PEAK, seqnames(R1PEAK))
R2PEAK <- split(R2PEAK, seqnames(R2PEAK))
R1PEAK <- R1PEAK[lengths(R1PEAK)>0]
R2PEAK <- R2PEAK[lengths(R2PEAK)>0]
chromosomes <- intersect(names(R1PEAK), names(R2PEAK))

countByOverlaps <- function(gi, pairs){
    counts_tab <- 0
    counts_r2 <- 0
    for(ps in pairs){
        ps <- read.delim(ps, comment.char="#", header=FALSE)
        ps <- with(ps, GInteractions(GRanges(V2, IRanges(V3, width=100)),
                                    GRanges(V4, IRanges(V5, width=100))))
        counts_tab <- counts_tab + countOverlaps(gi, ps, use.region="same")
        counts_r2 <- counts_r2 + countOverlaps(second(gi), second(ps))
    }
    list(count=counts_tab, shortCount=counts_r2)
}

getMscore <- function(mscore, gr){
    vw <- Views(mscore, unique(gr))
    sc <- viewMeans(vw)
    vs <- ranges(vw)
    vs <- as(vs, "GRanges")
    vs$score <- unlist(sc, use.names=FALSE)
    ol <- findOverlaps(gr, vs, type="equal")
    ol <- ol[!duplicated(queryHits(ol))]
    score <- vs[subjectHits(ol)]$score[match(seq_along(gr), queryHits(ol))]
    score[is.na(score)] <- 0
    score
}

### load counts
gis <- NULL
for(chrom in chromosomes){
    r1peak <- R1PEAK[[chrom]]
    r2peak <- R2PEAK[[chrom]]
    peak_pair <- expand.grid(seq_along(r1peak), seq_along(r2peak))
    gi <- GInteractions(r1peak[peak_pair[, 1]], r2peak[peak_pair[, 2]])
    cnt <- countByOverlaps(gi, pairs)
    gi$count <- cnt$count
    gi$shortCount <- cnt$shortCount
    gi <- gi[gi$count>0]
    if(length(gi)>0){
        gis <- c(gis, gi)
    }
}
gis <- do.call(c, gis)
### load gc content
fa <- FaFile(file=FASTA)
gc1 <- letterFrequency(getSeq(fa, first(gis)), letters="CG", as.prob=TRUE)
gc2 <- letterFrequency(getSeq(fa, second(gis)), letters="CG", as.prob=TRUE)
gis$gc <- gc1 * gc2 + 1e-9
### load enzyme cut number in R1
cut <- import(CUT)
start(cut) <- end(cut)
gis$cut <- countOverlaps(first(gis), cut)+0.1
### load mapping score
mscore <- import(MAPPABILITY, as="RleList")
m1 <- getMscore(mscore, first(gis))
m2 <- getMscore(mscore, second(gis))
gis$mappability <- m1*m2 + 1e-6
### get distance of the anchors
gis$dist <- distance(first(gis), second(gis))+1

mm <- log(as.data.frame(mcols(gis)[, c("cut", "gc", "mappability", "dist", "shortCount")]))
mm <- cbind(as.data.frame(first(gis)), as.data.frame(second(gis)), gis$count, mm)
colnames(mm) <- c("chr1", "start1", "end1", "width1", "strand1",
                "chr2", "start2", "end2", 'width2', "strand2",
                "count", "logl", "loggc", "logm", "logdist", 'logShortCount')

## doing statistics and resampling

pospoisson_regression <- function(mm) {
    fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
    # fit <- vglm(count ~ loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
    mm$expected = fitted(fit)
    mm$p_val = ppois(mm$count, mm$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected, lower.tail = FALSE, log.p = FALSE)
    m1 = mm[ mm$p_val > 1/length(mm$p_val),]
    fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
    # fit <- vglm(count ~  loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
    coeff<-round(coef(fit),10)
    mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10)
    # mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10)
    mm$expected2 <- mm$expected2 /(1-exp(-mm$expected2))
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = ppois(mm$count, mm$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected2, lower.tail = FALSE, log.p = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * nrow(mm)
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

negbinom_regression <- function(mm) {
    fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data = mm)
    # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = mm)
    mm$expected = fitted(fit)
    sze = fit$theta ##size parameter
    mm$p_val = pnbinom(mm$count, mu = mm$expected, size = sze, lower.tail = FALSE)
    m1 = mm[ mm$p_val > ( 1 / length(mm$p_val)),]
    ## second regression
    fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data = m1)
    # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = m1)
    coeff<-round(fit$coefficients,10)
    sze = fit$theta
    mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10) ## mu parameter
    # mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10) ## mu parameter
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = pnbinom(mm$count, mu = mm$expected2, size = sze, lower.tail = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * nrow(mm)
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
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
    group$weight <- 1
    group <- graphBAM(group)
    group <- connectedComp(ugraph(group))
    group <- lapply(group, as.numeric)
    group <- data.frame(id=unlist(group), g=rep(seq_along(group), lengths(group)))

    final$Cluster <- NA
    final$Cluster[group$id] <- group$g
    final$ClusterSize <- 0
    final$ClusterSize[group$id] <- table(group$g)[group$g]
    final$Cluster[is.na(final$Cluster)] <- seq(from=max(group$g)+1, to=max(group$g)+sum(is.na(final$Cluster)))
    final$NegLog10P <- -log10( final$p_val_reg2 )
    NegLog10P <- rowsum(final$NegLog10P, final$Cluster)
    final$NegLog10P <- NegLog10P[final$Cluster, 1]

    x <- unique( final[ final$ClusterSize != 0, c('chr1', 'Cluster', 'NegLog10P', 'ClusterSize')] )
    if(nrow(x)==0){
        final$ClusterType <- 'Singleton'
        return(final)
    }

    # sort rows by cumulative -log10 P-value
    x <- x[ order(x$NegLog10P) ,]
    y<-sort(x$NegLog10P)
    z<-cbind( seq(1,length(y),1), y )

    # keep a record of z before normalization
    z0 <- z

    z[,1]<-z[,1]/max(z[,1])
    z[,2]<-z[,2]/max(z[,2])

    u<-z
    u[,1] <-  1/sqrt(2)*z[,1] + 1/sqrt(2)*z[,2]
    u[,2] <- -1/sqrt(2)*z[,1] + 1/sqrt(2)*z[,2]

    v<-cbind(u, seq(1,nrow(u),1) )
    RefPoint <- v[ v[,2]==min(v[,2]) , 3]
    RefValue <- z0[RefPoint,2]

    # define peak cluster type
    final$ClusterType <- rep(NA, nrow(final))
    if(length(ol_)) final$ClusterType[ ol_ ] <- 'Singleton'
    if(length(ol)){
        final$ClusterType[ seq_along(gi) %in% ol & final$NegLog10P<RefValue  ] <-  'SharpPeak'
        final$ClusterType[ seq_along(gi) %in% ol & final$NegLog10P>=RefValue  ] <- 'BroadPeak'
    }
    return(final)
}

mm <- switch(REG_TYPE,
            "pospoisson"=pospoisson_regression(mm),
            "negbinom"=negbinom_regression(mm))

mm = classify_peaks(mm)

peaks <- if(nrow(mm)>0) subset(mm, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > FDR) else data.frame()
if (dim(peaks)[1] == 0) {
    print(paste('ERROR MAPS_regression_and_peak_caller.r: 0 bin pairs with count >= ',COUNT_CUTOFF,' observed/expected ratio >= ',RATIO_CUTOFF,' and -log10(fdr) > ',fdr_cutoff,sep=''))
    quit()
}

outf_name = paste(OUTPUT, '.',FDR,'.peaks',sep='')
dir.create(OUTPUT, recursive=TRUE)
write.table(peaks, file.path(OUTPUT, outf_name),
            row.names = FALSE, col.names = TRUE, quote=FALSE)

summary_all_runs <- split(peaks, peaks$ClusterType)
summary_all_runs <- lapply(summary_all_runs, function(.ele){
    c(count = nrow(.ele),
    minWidth1 = min(.ele$width1),
    medianWidth1 = median(.ele$width1),
    maxWidth1 = max(.ele$width1),
    minWidth2 = min(.ele$width2),
    medianWidth2 = median(.ele$width2),
    maxWidth2 = max(.ele$width2),
    minFoldChange = min(.ele$ratio2),
    medianFoldChange = median(.ele$ratio2),
    maxFoldChange = max(.ele$ratio2))
})
summary_all_runs <- do.call(rbind, summary_all_runs)
summary_outf_name = paste('summary.',OUTPUT,'.txt',sep='')
write.table(summary_all_runs, file.path(OUTPUT, summary_outf_name), row.names = TRUE, col.names = TRUE, quote=FALSE)
