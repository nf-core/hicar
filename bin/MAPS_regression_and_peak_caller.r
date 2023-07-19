#!/usr/bin/env Rscript

#########################################
# Author: [Ivan Juric](https://github.com/ijuric)
# File: MAPS_regression_and_peak_caller.r
# Source: https://github.com/ijuric/MAPS/blob/master/bin/MAPS/MAPS_regression_and_peak_caller.r
# Source+commit: https://github.com/ijuric/MAPS/blob/e6d5fdee7241b8f9a466cc778be6d0769b984e81/bin/MAPS/MAPS_regression_and_peak_caller.r
# Data: 11/08/2021, commit: e6d5fde
# This source code is licensed under the GPL-3.0 license
# modified by Jianhong:
# ## 1. set the Rscript environment.
# ## 2. export the COUNT_CUTOFF, RATIO_CUTOFF and FDR parameters
# ## 3. Automatic detect the chromosome names
# ## 4. Handle the error if the input count table is empty
# ## 5. Prefilter the data before fit to vglm to handle the NA error by replace the loglikelihood function.
# ## 6. Handle the error if output is empty
# ## 7. Handle the error if AND or XOR table is empty
# ## 8. clean unused code
# ## 9. fix the indent space
#########################################

## run example:
##  Rscript MAPS_regression_and_peak_caller.r /home/jurici/work/PLACseq/MAPS_pipe/results/mESC_test/ MY_115.5k 5000 1 None pospoisson NA
##
## arguments:
## INFDIR - dir with reg files
## SET - dataset name
## RESOLUTION - resolution (for example 5000 or 10000)
## COUNT_CUTOFF - count cutoff, default 12
## RATIO_CUTOFF - ratio cutoff, default 2.0
## FDR - -log10(fdr) cutoff, default 2
## FILTER - file containing bins that need to be filtered out. Format: two columns "chrom", "bin". "chrom" contains 'chr1','chr2',.. "bin" is bin label
## regresison_type - pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson

library(VGAM)
library(MASS)
options(warn=-1)

### constants
chroms = NULL
runs = c(1)
RESOLUTION = NULL

COUNT_CUTOFF = 12
RATIO_CUTOFF = 2.0
GAP = 15000
FDR = 2

REG_TYPE = 'pospoisson'
###

args <- commandArgs(trailingOnly=TRUE)
fltr = data.frame(chr='chrNONE',bin=-1)

if (length(args) < 3 || length(args) > 8) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): INFDIR, SET, RESOLUTION, COUNT_CUTOFF, RATIO_CUTOFF, FDR, FILTER, regression_type.')
    print('FILTER is optional argument. Omitt it if no filtering required.')
    print(paste('Number of arguments entered:',length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    INFDIR = args[1]
    SET = args[2]
    RESOLUTION = as.integer(args[3])
    SET <- paste0(SET, ".", ceiling(RESOLUTION/1e3), "k")
    chroms <- dir(INFDIR, "reg_raw.*")
    chroms <- unique(sub("reg_raw\\.(.*?)\\..*$", "\\1", chroms))
    chroms <- unique(unlist(strsplit(chroms, "_")))
    if(length(args)>3){
        COUNT_CUTOFF = as.numeric(args[4])
        RATIO_CUTOFF = as.numeric(args[5])
        FDR = as.numeric(args[6])
        print('filter used (if any):')
        if (length(args) > 6) {
            if (args[7] != 'None') {
                FILTER = args[7]
                fltr = read.table(FILTER,header=T)
                print(fltr)
            } else {
                print('None')
            }
            ## this done so that the script is compatible with previous run_pipeline scripts
            if (length(args) > 7) {
                if (args[8] != 'pospoisson' && args[8] != 'negbinom') {
                    print(paste('wrong regression choice. Your choice:', args[8], '. Avaiable choices: pospoisson or negbinom'),sep = ' ')
                    quit()
                }
                REG_TYPE = args[8]
            }
        } else {
            print('None')
        }
    }
}

## loading data
mm_combined_and = data.frame()
mm_combined_xor = data.frame()
outf_names = c()
for (i in chroms) {
    for (k in chroms){
        for (j in c('.and','.xor')) {
            print(paste('loading chromosome ',i,' ', k, ' ',j,sep=''))
            inf_name = paste(INFDIR,'reg_raw.',i,'_', k, '.',SET,sep='')
            if(file.exists(paste(inf_name,j,sep=''))){
                outf_names = c(outf_names, paste(inf_name,j,'.MAPS2_',REG_TYPE,sep = ''))
                mm = read.table(paste(inf_name,j,sep=''),header=T)
                mm$chr = rep(i, nrow(mm))
                mm$chr2 = rep(k, nrow(mm))
                mm = subset( mm, dist > 1) # removing adjacent bins
                mm = subset(mm, !(mm$chr %in% fltr$chr & (mm$bin1_mid %in% fltr$bin | mm$bin2_mid %in% fltr$bin ))) ## filtering out bad bins
                if (j == '.and') {
                    mm_combined_and = rbind(mm_combined_and, mm)
                } else if (j == '.xor') {
                    mm_combined_xor = rbind(mm_combined_xor, mm)
                }
            }
        }
    }
}

dataset_length_and = length(mm_combined_and$bin1_mid)
dataset_length_xor = length(mm_combined_xor$bin1_mid)
dataset_length = dataset_length_and + dataset_length_xor

## doing statistics and resampling
loglikelihood <- function (mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
    lambda <- eta2theta(eta, "loglink", earg = list(bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE,
        tag = FALSE))
    if (residuals) {
        stop("loglikelihood residuals not implemented yet")
    }
    else {
        if(exists('dgaitdpois', where = 'package:VGAM', mode='function')){
            dgaitpois <- dgaitdpois
        }
        ll.elts <- c(w) * dgaitpois(y, lambda, truncate = 0,
            log = TRUE)
        if (summation) {
            sum(ll.elts[!is.infinite(ll.elts)])
        }
        else {
            ll.elts
        }
    }
}
pospoisson_regression <- function(mm, dataset_length) {
    family <- pospoisson()
    family@loglikelihood <- loglikelihood
    # fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
    fit <- vglm(count ~ loggc + logm + logdist + logShortCount, family = family, data = mm)
    mm$expected = fitted(fit)
    mm$p_val = ppois(mm$count, mm$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected, lower.tail = FALSE, log.p = FALSE)
    m1 = mm[ mm$p_val > 1/length(mm$p_val),]
    # fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
    fit <- vglm(count ~  loggc + logm + logdist + logShortCount, family = family, data = m1)
    coeff<-round(coef(fit),10)
    # mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10)
    mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10)
    mm$expected2 <- mm$expected2 /(1-exp(-mm$expected2))
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = ppois(mm$count, mm$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected2, lower.tail = FALSE, log.p = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * dataset_length
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

negbinom_regression <- function(mm, dataset_length) {
    #fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data = mm)
    fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = mm)
    mm$expected = fitted(fit)
    sze = fit$theta ##size parameter
    mm$p_val = pnbinom(mm$count, mu = mm$expected, size = sze, lower.tail = FALSE)
    m1 = mm[ mm$p_val > ( 1 / length(mm$p_val)),]
    ## second regression
    # fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data = m1)
    fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = m1)
    coeff<-round(fit$coefficients,10)
    sze = fit$theta
    #mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10) ## mu parameter
    mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10) ## mu parameter
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = pnbinom(mm$count, mu = mm$expected2, size = sze, lower.tail = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * dataset_length
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

do_summaries <- function(peaks_and,peaks_xor, peaks, fraction, r) {
    ## no peaks with this fdr
    if (ncol(peaks_and) == 0) {
        peaks_and = data.frame('count' = NA, 'dist'=NA, 'p_val_reg2'=NA, 'fdr'=NA)
    }
    if (ncol(peaks_xor) == 0) {
        peaks_xor = data.frame('count' = NA, 'dist'=NA, 'p_val_reg2'=NA, 'fdr'=NA)
    }
    if (ncol(peaks_and) == 0 & ncol(peaks_xor) == 0) {
        peaks = data.frame('count' = NA, 'dist'=NA, 'p_val_reg2'=NA, 'fdr'=NA)
    }
    summary_one_fdr_val = data.frame('run' = r,
                                    'log10_fdr_cutoff' = fdr_cutoff,
                                    'count_cutoff' = COUNT_CUTOFF,
                                    'ratio_cutoff' = RATIO_CUTOFF,
                                    'total_loops' = nrow(peaks),
                                    'total_AND_loops'=nrow(peaks_and),
                                    'total_XOR_loops'=nrow(peaks_xor),
                                    'singleton_fraction' = fraction,
                                    'AND_size'=length(peaks_and$count),
                                    'AND_mean_dist'=mean(peaks_and$dist, na.rm=TRUE)*RESOLUTION,
                                    'AND_median_dist'=median(peaks_and$dist, na.rm=TRUE)*RESOLUTION,
                                    'AND_min_count'=min(peaks_and$count, na.rm=TRUE),
                                    'AND_max_pval'=max(peaks_and$p_val_reg2, na.rm=TRUE),
                                    'AND_max_fdr'=max(peaks_and$fdr, na.rm=TRUE),
                                    'XOR_size'=length(peaks_xor$count),
                                    'XOR_mean_dist'=mean(peaks_xor$dist, na.rm=TRUE)*RESOLUTION,
                                    'XOR_median_dist'=median(peaks_xor$dist, na.rm=TRUE)*RESOLUTION,
                                    'XOR_min_count'=min(peaks_xor$count, na.rm=TRUE),
                                    'XOR_max_pval'=max(peaks_xor$p_val_reg2, na.rm=TRUE),
                                    'XOR_max_fdr'=max(peaks_xor$fdr, na.rm=TRUE),
                                    'size'=length(peaks$count),
                                    'mean_dist'=mean(peaks$dist, na.rm=TRUE)*RESOLUTION,
                                    'median_dist'=median(peaks$dist, na.rm=TRUE)*RESOLUTION,
                                    'min_count'=min(peaks$count, na.rm=TRUE),
                                    'max_pval'=max(peaks$p_val_reg2, na.rm=TRUE),
                                    'max_fdr'=max(peaks$fdr, na.rm=TRUE))
    return(summary_one_fdr_val)
}

label_peaks <- function(df) {
    chroms = unique(df$chr)
    print('chromosomes with potential interactions:')
    print(chroms)
    final = data.frame()
    for (CHR in chroms) {
        y = df[df$chr == CHR,]
        y$p_val_reg2[ y$p_val_reg2 == 0 ] = 1111111
        y$p_val_reg2[ y$p_val_reg2 == 1111111 ] = min(y$p_val_reg2)
        for(i in 1:nrow(y)) {
            z <- y[ abs(y$bin1_mid - y$bin1_mid[i])<=GAP & abs(y$bin2_mid - y$bin2_mid[i])<=GAP,]
            y$CountNei[i] <- nrow(z)
        }
        u <- y[ y$CountNei == 1 ,] # singletons
        v <- y[ y$CountNei >= 2 ,] # peak cluster: sharp peak + broad peak
        out <- NULL
        if(nrow(u)>0) {
            u$label <- 0 # for singletons, assign cluster label 0
            # for singletons, cluster size = 1
            u$NegLog10P <- -log10(u$p_val_reg2 )
            u$ClusterSize <- 1
            out <- rbind(out, u)
        }
        if(nrow(v)>0) {
            v$label <- seq(1,nrow(v),1) # for peak cluster, assign label 1, 2, ..., N
            # for all bin pairs within the neighborhood
            # assign the same cluster label, using the minimal label
            for(i in 1:nrow(v)) {
                w <- v[ abs(v$bin1_mid - v$bin1_mid[i])<=GAP & abs(v$bin2_mid - v$bin2_mid[i])<=GAP ,]
                w.min <- min(w$label)
                w.label <- sort(unique(w$label))
                for(j in 2:length(w.label)) {
                    v$label[ v$label == w.label[j] ] <- w.min
                }
            }
            # step 4: assign consecutive label number
            v.rec <- sort( unique( v$label ) )
            v.rec <- cbind(v.rec, seq(1, length(v.rec), 1))
            for(i in 1:nrow(v)) {
                v$label[i] <- v.rec[ v.rec[,1]==v$label[i] ,2]
            }
            # step 5: calculate cumulative NegLog10P and the cluster size
            v$NegLog10P <- 0
            v$ClusterSize <- 0
            for(i in 1:nrow(v.rec)) {
            # find all bin pairs within the same cluster i
                vtmp <- v[ v$label == i ,]
                v$NegLog10P[ v$label == i ] <- sum( -log10( vtmp$p_val_reg2 ) )
                v$ClusterSize[ v$label == i ] <- nrow(vtmp)
            }
            out <- rbind(out, v)
        }
        final<-rbind(final, out)
    }
    print(dim(final))
    return(final)
}

classify_peaks <- function(final) {

    # only keep unique peak clusters, not bin pairs
    x <- unique( final[ final$label != 0, c('chr', 'label', 'NegLog10P', 'ClusterSize')] )
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
    RefPoint <- v[ v[,2]==min(v[,2]) , 3] # 1423
    RefValue <- z0[RefPoint,2]

    # define peak cluster type
    final$ClusterType <- '0'
    final$ClusterType[ final$label==0 ] <- 'Singleton'
    final$ClusterType[ final$label>=1 & final$NegLog10P<RefValue  ] <-  'SharpPeak'
    final$ClusterType[ final$label>=1 & final$NegLog10P>=RefValue  ] <- 'BroadPeak'
    #table(final$ClusterType)
    return(final)
}


mx_combined_and = data.frame()
mx_combined_xor = data.frame()
summary_all_runs = data.frame()

singletons_names = paste(chroms,'_0',sep='')
for (r in runs) {
## in case you want to do resampling
## you'd put resampling code here
    for (i in chroms) {
        for(j in chroms) {
            ## regression
            tryCatch({
                print(paste('run',r,': regression on chromosome',i, ' and ', j))
                outf_names <- paste(INFDIR, 'reg_raw.',i,'_', j, '.and.',SET,'.MAPS2_',REG_TYPE,sep = '')
                mm = subset(mm_combined_and, chr == i & chr2 == j)
                if(nrow(mm)>6){
                    if (REG_TYPE == 'pospoisson') {
                        mm = pospoisson_regression(mm, dataset_length)
                    } else if (REG_TYPE == 'negbinom') {
                        mm = negbinom_regression(mm, dataset_length)
                    }
                    mx_combined_and = rbind(mx_combined_and, mm)
                }
                write.table(mm,outf_names,row.names = TRUE,col.names = TRUE,quote=FALSE)
                outf_names <- paste(INFDIR, 'reg_raw.',i,'_', j, '.xor.',SET,'.MAPS2_',REG_TYPE,sep = '')
                mm = subset(mm_combined_xor, chr == i & chr2 == j)
                if(nrow(mm)>6){
                    if (REG_TYPE == 'pospoisson') {
                        mm = pospoisson_regression(mm, dataset_length)
                    } else if (REG_TYPE == 'negbinom') {
                        mm = negbinom_regression(mm, dataset_length)
                    }
                    mx_combined_xor = rbind(mx_combined_xor, mm)
                }
                write.table(mm,outf_names,row.names = TRUE,col.names = TRUE,quote=FALSE)
            }, error=function(.e){
                message(.e)
            })
            ## finding min FDR so I can set appropriate lower boundary for FDR
        }
    }
    ## save QC file
    qc_out = paste(INFDIR,SET,'.maps.qc',sep = '')
    qc_label = c('AND_set','XOR_set')
    qc_val = c(sum(mm_combined_and$count), sum(mm_combined_xor$count))
    qc_name = c('number of sequencing pairs in AND set', 'number of sequencing pairs in XOR set')
    df_qc = data.frame()

    summary_one_run = data.frame()
    singletons = data.frame()
    for (fdr_cutoff in FDR) {
        peaks_and = if(nrow(mx_combined_and)>0) subset(mx_combined_and, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > fdr_cutoff) else data.frame()
        peaks_xor = if(nrow(mx_combined_xor)>0) subset(mx_combined_xor, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > fdr_cutoff) else data.frame()
        print("finding peaks")
        peaks = rbind(peaks_and, peaks_xor)
        if (dim(peaks)[1] == 0) {
            print(paste('ERROR MAPS_regression_and_peak_caller.r: 0 bin pairs with count >= ',COUNT_CUTOFF,' observed/expected ratio >= ',RATIO_CUTOFF,' and -log10(fdr) > ',fdr_cutoff,sep=''))
            quit()
        }
        peaks = label_peaks(peaks)
        peaks$lab = paste(peaks$chr, peaks$label,sep='_')
        peak_types = classify_peaks(peaks)
        outf_name = paste(INFDIR,SET, '.',fdr_cutoff,'.peaks',sep='')
        write.table(peak_types,outf_name, row.names = FALSE, col.names = TRUE, quote=FALSE)
        print("finding singletons")
        peak_classes = table(peaks$lab)
        n_singletons = sum(peak_classes[names(peak_classes) %in% singletons_names])
        n_singleton_chroms = length(peak_classes[names(peak_classes) %in% singletons_names])

        fraction= n_singletons / (length(peak_classes) + n_singletons - n_singleton_chroms)
        fraction = n_singletons / (length(peak_classes) + n_singletons - n_singleton_chroms)

        singletons_one_run = data.frame('fdr'=fdr_cutoff, 'fraction'=fraction)
        singletons = rbind(singletons, singletons_one_run)
        print(paste(fdr_cutoff,':',singletons))
        summary_one_run = rbind(summary_one_run, do_summaries(peaks_and, peaks_xor, peaks, fraction, r))
    }
    ## find singletons, sharp peaks, broad peaks
    summary_all_runs = rbind(summary_all_runs, summary_one_run)
}
summary_outf_name = paste(INFDIR,'summary.',SET,'.txt',sep='')
write.table(summary_all_runs, summary_outf_name, row.names = FALSE, col.names = TRUE, quote=FALSE)
