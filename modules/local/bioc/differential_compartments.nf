process DIFFERENTIAL_COMPARTMENTS {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::bioconductor-diffhic=1.26.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-diffhic:1.26.0--r41hc247a5b_2' :
        'biocontainers/bioconductor-diffhic:1.26.0--r41hc247a5b_2' }"

    input:
    tuple val(bin_size), val(meta), path(bigwigs, stageAs: "adjust_compartments/*")

    output:
    tuple val(bin_size), path("${prefix}_pval.bigWig")                                     , emit: diff
    path "versions.yml"                                                                    , emit: versions

    script:
    prefix   = task.ext.prefix ?: "diffCompartments_bin${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Feb. 22, 2024 do differential analysis for the A/B compartments
    ## the compartment score were fist normalized by limma::normalizeQuantiles
    ## P-Value is calculated by chiq-square test for Mahalanobis Distance
    ## Copyright (c) 2024 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("limma", "rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    OUTFILE = paste0("$prefix", "_pval.bigWig")
    fs <- dir('adjust_compartments', '.(bigWig|bw)', full.names=TRUE, ignore.case = TRUE)
    names(fs) <- basename(fs)

    ## import all A/B compartment scores
    d <- lapply(fs, import)
    if(length(fs)<2){
        out <- lapply(d, function(.ele){
            mcols(.ele)[, 'score'] <- 1
            .ele
        })
    }else{
        ## split the regions into non-overlapping regions
        regions <- disjoin(unlist(GRangesList(d)))
        ol <- lapply(d, function(.ele){
            .ol <- findOverlaps(regions, .ele, minoverlap=1L)
            .id <- subjectHits(.ol)[match(seq_along(regions), queryHits(.ol))]
            out <- rep(0, length(regions))
            out[!is.na(.id)] <- mcols(.ele[.id[!is.na(.id)]])[, 'score']
            out
        })
        mcols(regions) <- do.call(cbind, ol)
        ## calculate p-value for each chromosome
        regions <- split(regions, seqnames(regions))
        dfs <- lapply(regions, function(.ele) as.data.frame(mcols(.ele)))
        ## normalization
        dfs <- lapply(dfs, function(df){
            df <- limma::normalizeQuantiles(df, ties=TRUE)
            df <- df/colMeans(abs(df), na.rm=TRUE)
        })
        ## calculate the centers
        intra_cen <- lapply(dfs, colMeans)
        ## cov
        intra_cov <- lapply(dfs, cov)
        ## number of samples
        n <- length(d)
        ## Mahalanobis distance
        maha <- mapply(dfs, intra_cen, intra_cov, FUN = function(df, cen, cov){
            apply(df, 1, function(.ele){
                .dif <- abs(diff(range(.ele)))^2
                if(.dif==0) return(NA)
                .cov <- solve(cov) %*% diag(.dif, n)
                mahalanobis(t(.ele), as.matrix(cen), .cov)
            })
        }, SIMPLIFY = FALSE)
        pval <- lapply(maha, function(.ele){
            p <- pchisq(.ele, df = n - 1, lower.tail = TRUE)
            p[is.na(p)] <- 1
            p
        })
        out <- mapply(regions, pval, FUN=function(gr, p){
            mcols(gr) <- data.frame(pval=p)
            gr
        })
        out <- unlist(GRangesList(out))
    }
    #mcols(out)[, 'score'] <- p.adjust(mcols(out)[, 'pval'], method = 'BH')
    mcols(out)[, 'score'] <- -10*log10(mcols(out)[, 'pval'])
    ## export pvalue for the difference of A/B compartment
    export(out[!is.na(mcols(out)[, 'score'])], OUTFILE)
    regions <- unlist(GRangesList(regions))
    ol <- findOverlaps(regions, out, type='equal')
    mcols(regions)[, "pval"] <- 1
    mcols(regions)[queryHits(ol), "pval"] <- mcols(out)[subjectHits(ol), 'pval']
    mcols(regions)[, 'score'] <- -10*log10(mcols(regions)[, 'pval'])
    rle <- list()
    for(n in colnames(mcols(regions))[seq_along(d)]){
        mcols(regions)[, paste0(n, '_AB')] <- ifelse(mcols(regions)[, n] >= 0, 'A', 'B')
        rle[[paste0(n, '_AB')]] <- mcols(regions)[, paste0(n, '_AB')]
    }
    write.csv(regions, sub('(bigWig|bw)', 'merged.csv', OUTFILE, ignore.case = TRUE))
    tab <- table(mcols(regions)[paste0(colnames(mcols(regions))[seq_along(d)], '_AB')])
    data.frame(tab)
    write.csv(data.frame(tab), sub('(bigWig|bw)', 'AB.tab.csv', OUTFILE, ignore.case = TRUE))
    if(length(rle)==2){
        rle <- do.call(cbind, rle)
        prle <- apply(rle, 1, paste, collapse='')
        prle <- rle(prle)
        prle <- table(prle[['values']])
        prle <- data.frame(do.call(rbind, strsplit(names(prle), split = '')),
                            prle)
        colnames(prle)[c(1, 2)] <- colnames(rle)
        prle <- prle[, c(colnames(rle), 'Freq')]
        write.csv(prle, sub('(bigWig|bw)', 'AB.rle.tab.csv', OUTFILE))
    }
    """
}
