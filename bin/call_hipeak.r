#!/usr/bin/env Rscript
#######################################################################
#######################################################################
## Created on Aug. 24, 2021 call peaks
## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
library(VGAM)
library(MASS)
writeLines(as.character(packageVersion("VGAM")), "VGAM.version.txt")
writeLines(as.character(packageVersion("MASS")), "MASS.version.txt")

options(warn=-1)

### constants
OUTPUT = "peaks.csv"
REG_TYPE = 'pospoisson'

if("optparse" %in% installed.packages()){
    library(optparse)
    option_list <- list(make_option(c("-m", "--regression_type"), type="character", default='pospoisson', help="pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson", metavar="string"),
                        make_option(c("-t", "--count"), type="character", default=NULL, help="count table output by prepare count", metavar="string"),
                        make_option(c("-o", "--output"), type="character", default="peaks.csv", help="output folder", metavar="string"))
    opt_parser <- OptionParser(option_list=option_list)
    opt <- parse_args(opt_parser)
}else{
    args <- commandArgs(TRUE)
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
    option_list <- list("regression_type"=c("--regression_type", "-m", "character"),
                        "count"=c("--count", "-t", "character"),
                        "output"=c("--output", "-o", "character"))
    opt <- parse_args(option_list, args)
}

if(!is.null(opt$output)){
    OUTPUT <- opt$output
}
if(!is.null(opt$output)){
    mm <- read.csv(opt$count)
}else{
    stop("count is required")
}

if(!all(c("chr1", "start1", "end1", "width1",
        "chr2", "start2", "end2", 'width2',
        "count", "logl", "logn", "loggc", "logm", "logdist", 'logShortCount') %in% colnames(mm))){
    stop("count table is not in correct format.")
}

## doing statistics and resampling

pospoisson_regression <- function(mm) {
    fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount + logn, family = pospoisson(), data = mm)
    # fit <- vglm(count ~ loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
    mm$expected = fitted(fit)
    mm$p_val = ppois(mm$count, mm$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected, lower.tail = FALSE, log.p = FALSE)
    m1 = mm[ mm$p_val > 1/length(mm$p_val),]
    fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount + logn, family = pospoisson(), data = m1)
    # fit <- vglm(count ~  loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
    coeff<-round(coef(fit),10)
    mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount + coeff[7]*mm$logn), 10)
    # mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10)
    mm$expected2 <- mm$expected2 /(1-exp(-mm$expected2))
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = ppois(mm$count, mm$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm$expected2, lower.tail = FALSE, log.p = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * nrow(mm)
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

negbinom_regression <- function(mm) {
    fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount + logn, data = mm)
    # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = mm)
    mm$expected = fitted(fit)
    sze = fit$theta ##size parameter
    mm$p_val = pnbinom(mm$count, mu = mm$expected, size = sze, lower.tail = FALSE)
    m1 = mm[ mm$p_val > ( 1 / length(mm$p_val)),]
    ## second regression
    fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount + logn, data = m1)
    # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = m1)
    coeff<-round(fit$coefficients,10)
    sze = fit$theta
    mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm + coeff[5]*mm$logdist + coeff[6]*mm$logShortCount) + coeff[7]*mm$logn, 10) ## mu parameter
    # mm$expected2 <- round(exp(coeff[1]  + coeff[2]*mm$loggc + coeff[3]*mm$logm + coeff[4]*mm$logdist + coeff[5]*mm$logShortCount), 10) ## mu parameter
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = pnbinom(mm$count, mu = mm$expected2, size = sze, lower.tail = FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * nrow(mm)
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

mm <- switch(REG_TYPE,
            "pospoisson"=pospoisson_regression(mm),
            "negbinom"=negbinom_regression(mm))

write.csv(mm, OUTPUT, row.names=FALSE)
