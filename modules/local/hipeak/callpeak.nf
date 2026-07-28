process CALL_HIPEAK {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bioconductor-monocle=2.20.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-monocle:2.20.0--r41h399db7b_0' :
        'biocontainers/bioconductor-monocle:2.20.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}.peaks.csv"), emit: peak
    path "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 call peaks
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    #######################################################################
    #######################################################################
    pkgs <- c("VGAM", "MASS")
    versions <- c("${task.process}:")
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
    if(grepl("-m", "$args") || grepl("--regression_type", "$args")){
        args <- strsplit("$args", "\\\\s+")[[1]]
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
    getFormula <- function(mm){
        coln <- c("logl", "loggc", "logm", "logdist", "logShortCount", "logn")
        ln <- vapply(coln, FUN=function(.ele){
            length(unique(mm[, .ele]))>1
        }, FUN.VALUE=logical(1))
        list("formula"=paste("count ~", paste(coln[ln], collapse="+")),
            "factor"=coln[ln])
    }
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
    pospoisson_regression <- function(mm) {
        dataset_length<- nrow(mm)
        formula <- getFormula(mm)
        fac <- formula[["factor"]]
        formula <- formula[["formula"]]
        family <- pospoisson()
        family@loglikelihood <- loglikelihood
        fit <- vglm(formula=as.formula(formula), family = family, data = mm)
        mm\$expected = fitted(fit)
        mm\$p_val = ppois(mm\$count, mm\$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm\$expected, lower.tail = FALSE, log.p = FALSE)
        m1 = mm[ mm\$p_val > 1/length(mm\$p_val),]
        fit <- vglm(formula=as.formula(formula), family = family, data = m1)
        coeff<-round(coef(fit),10)
        mm\$expected2 <- round(exp(coeff[1] + rowSums(t(coeff[-1]*t(mm[, fac])))), 10)
        mm\$expected2 <- mm\$expected2 /(1-exp(-mm\$expected2))
        mm\$ratio2 <- mm\$count / mm\$expected2
        mm\$p_val_reg2 = ppois(mm\$count, mm\$expected2, lower.tail = FALSE, log.p = FALSE) / ppois(0, mm\$expected2, lower.tail = FALSE, log.p = FALSE)
        mm\$p_bonferroni = mm\$p_val_reg2 * dataset_length
        mm\$fdr <- p.adjust(mm\$p_val_reg2, method='fdr')
        return(mm)
    }

    negbinom_regression <- function(mm) {
        formula <- getFormula(mm)
        fac <- formula[["factor"]]
        formula <- formula[["formula"]]
        fit <- glm.nb(formula=as.formula(formula), data = mm)
        mm\$expected = fitted(fit)
        sze = fit\$theta ##size parameter
        mm\$p_val = pnbinom(mm\$count, mu = mm\$expected, size = sze, lower.tail = FALSE)
        m1 = mm[ mm\$p_val > ( 1 / length(mm\$p_val)),]
        ## second regression
        fit <- glm.nb(formula=as.formula(formula), data = m1)
        coeff<-round(fit\$coefficients,10)
        sze = fit\$theta
        mm\$expected2 <- round(exp(coeff[1] + rowSums(t(coeff[-1]*t(mm[, fac])))), 10) ## mu parameter
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
