// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CALL_HIPEAK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-monocle=2.20.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-monocle:2.20.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-monocle:2.20.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.peaks.csv"), emit: peak
    path "versions.yml"                          , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 call peaks
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("VGAM", "MASS")
    versions <- c("${getProcessName(task.process)}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # wirte versions.yml

    options(warn=-1)
    ### Options
    ## make_option(c("-m", "--regression_type"), type="character", default='pospoisson', help="pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson", metavar="string"),
    ## make_option(c("-t", "--count"), type="character", default=NULL, help="count table output by prepare count", metavar="string"),
    ## make_option(c("-o", "--output"), type="character", default="peaks.csv", help="output folder", metavar="string")

    OUTPUT = "${meta.id}.peaks.csv"
    REG_TYPE = 'pospoisson'
    mm <- read.csv("$counts")
    if(grepl("-m", "$options.args") || grepl("--regression_type", "$options.args")){
        args <- strsplit("$options.args", "\\\\s+")[[1]]
        id <- args=="-m" | args=="--regression_type"
        id <- which(id)[1]
        if(args[id+1] %in% c("negbinom", "pospoisson")){
            REG_TYPE <- args[id+1]
        }
    }

    if(!all(c("chr1", "start1", "end1", "width1",
            "chr2", "start2", "end2", 'width2',
            "count", "logl", "logn", "loggc", "logm", "logdist", 'logShortCount') %in% colnames(mm))){
        stop("count table is not in correct format.")
    }

    ## doing statistics and resampling
    checkdata <- function(mm){
        mf <- vglm(count ~ logl + loggc + logm + logdist + logShortCount + logn, family = pospoisson(), data = mm, method="model.frame")
        y <- model.response(mf, "any")
        w <- model.weights(mf)
        if (!length(w)) {
            w <- rep_len(1, nrow(mf))
        }
        lambda.init <- Init.mu(y = y, w = w, imethod = 1, imu = NULL)
        eta <- theta2eta(lambda.init, "loglink", earg = list(
                theta = NULL, bvalue = NULL, inverse = FALSE, deriv = 0,
                short = TRUE, tag = FALSE))
        lambda <- eta2theta(eta, "loglink", earg = list(theta = NULL,
            bvalue = NULL, inverse = FALSE, deriv = 0, short = TRUE,
            tag = FALSE))
        ll.elts <- dgaitpois(y, lambda, truncate = 0, log = TRUE)
    }
    trimData <- function(mm){
        ll.elts <- checkdata(mm)
        while(any(is.infinite(ll.elts))){
            mm <- mm[!is.infinite(ll.elts), , drop=FALSE]
            ll.elts <- checkdata(mm)
        }
        mm
    }
    pospoisson_regression <- function(mm) {
        dataset_length<- nrow(mm)
        mm <- trimData(mm)
        fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount + logn, family = pospoisson(), data = mm)
        # fit <- vglm(count ~ loggc + logm + logdist + logShortCount, family = pospoisson(), data = mm)
        mm\$expected = fitted(fit)
        mm\$p_val = ppois(mm\$count, mm\$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm\$expected, lower.tail = FALSE, log.p = FALSE)
        m1 = mm[ mm\$p_val > 1/length(mm\$p_val),]
        m1 <- trimData(m1)
        fit <- vglm(count ~ logl + loggc + logm + logdist + logShortCount + logn, family = pospoisson(), data = m1)
        # fit <- vglm(count ~  loggc + logm + logdist + logShortCount, family = pospoisson(), data = m1)
        coeff<-round(coef(fit),10)
        mm\$expected2 <- round(exp(coeff[1] + coeff[2]*mm\$logl + coeff[3]*mm\$loggc + coeff[4]*mm\$logm + coeff[5]*mm\$logdist + coeff[6]*mm\$logShortCount + coeff[7]*mm\$logn), 10)
        # mm\$expected2 <- round(exp(coeff[1]  + coeff[2]*mm\$loggc + coeff[3]*mm\$logm + coeff[4]*mm\$logdist + coeff[5]*mm\$logShortCount), 10)
        mm\$expected2 <- mm\$expected2 /(1-exp(-mm\$expected2))
        mm\$ratio2 <- mm\$count / mm\$expected2
        mm\$p_val_reg2 = ppois(mm\$count, mm\$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm\$expected2, lower.tail = FALSE, log.p = FALSE)
        mm\$p_bonferroni = mm\$p_val_reg2 * dataset_length
        mm\$fdr <- p.adjust(mm\$p_val_reg2, method='fdr')
        return(mm)
    }

    negbinom_regression <- function(mm) {
        fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount + logn, data = mm)
        # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = mm)
        mm\$expected = fitted(fit)
        sze = fit\$theta ##size parameter
        mm\$p_val = pnbinom(mm\$count, mu = mm\$expected, size = sze, lower.tail = FALSE)
        m1 = mm[ mm\$p_val > ( 1 / length(mm\$p_val)),]
        ## second regression
        fit <- glm.nb(count ~ logl + loggc + logm + logdist + logShortCount + logn, data = m1)
        # fit <- glm.nb(count ~  loggc + logm + logdist + logShortCount, data = m1)
        coeff<-round(fit\$coefficients,10)
        sze = fit\$theta
        mm\$expected2 <- round(exp(coeff[1] + coeff[2]*mm\$logl + coeff[3]*mm\$loggc + coeff[4]*mm\$logm + coeff[5]*mm\$logdist + coeff[6]*mm\$logShortCount) + coeff[7]*mm\$logn, 10) ## mu parameter
        # mm\$expected2 <- round(exp(coeff[1]  + coeff[2]*mm\$loggc + coeff[3]*mm\$logm + coeff[4]*mm\$logdist + coeff[5]*mm\$logShortCount), 10) ## mu parameter
        mm\$ratio2 <- mm\$count / mm\$expected2
        mm\$p_val_reg2 = pnbinom(mm\$count, mu = mm\$expected2, size = sze, lower.tail = FALSE)
        mm\$p_bonferroni = mm\$p_val_reg2 * nrow(mm)
        mm\$fdr <- p.adjust(mm\$p_val_reg2, method='fdr')
        return(mm)
    }

    mm <- switch(REG_TYPE,
                "pospoisson"=pospoisson_regression(mm),
                "negbinom"=negbinom_regression(mm))

    write.csv(mm, OUTPUT, row.names=FALSE)
    """
}
