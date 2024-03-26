process ADJUST_COMPARTMENTS {
    tag "$bin_size"
    label 'process_medium'

    conda "bioconda::bioconductor-trackviewer=1.38.1"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.38.1--r43hdfd78af_0' :
        'biocontainers/bioconductor-trackviewer:1.38.1--r43hdfd78af_0' }"

    input:
    tuple val(bin_size), val(meta), path(bigwigs, stageAs: "old_compartments/*")

    output:
    tuple val(bin_size), val(meta), path("${prefix}/*")                                    , emit: compartments
    path "versions.yml"                                                                    , emit: versions

    script:
    prefix   = task.ext.prefix ?: "adjust_compartments_bin${bin_size}"
    """
    #!/usr/bin/env Rscript
    #######################################################################
    #######################################################################
    ## Created on Feb. 22, 2024 adjust the A/B compartments by control group
    ## The hypothesis for this adjustment is that most of the compartment
    ##  should keep same.
    ## Copyright (c) 2024 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("trackViewer", "rtracklayer")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    OUTFOLDER = "$prefix"
    dir.create(OUTFOLDER, showWarnings = FALSE)
    fs <- dir('old_compartments', '.(bigWig|bw)', full.names=TRUE, ignore.case = TRUE)
    names(fs) <- basename(fs)
    if(length(fs)<2){
        file.copy(fs, to = file.path(OUTFOLDER, basename(fs)))
    }else{
        ## import all A/B compartment scores
        d <- lapply(fs, import)
        ## split the regions into non-overlapping regions
        regions <- disjoin(unlist(GRangesList(d)))
        mcols(regions)[, "score"] <- 0
        ## get the signature for each sample for the disjoin regions
        regions_sign <- lapply(d, function(.ele){
            mcols(GRoperator(regions, .ele, operator = '+'))[, "score"]
        })
        regions_sign <- do.call(cbind, regions_sign)
        regions_sign <- sign(regions_sign)
        regions_sign[is.na(regions_sign)] <- 1
        ## get the product for each row, if over 50% sign are same, the product should be positive.
        ## eg -1 X -1 = 1; 1 X 1 = 1
        mcols(regions) <- regions_sign
        regions <- split(regions, as.character(seqnames(regions)))
        regions_sign_summary <- lapply(regions, function(.ele) {
            if(length(.ele)<1){
                return(FALSE)
            }
            mc <- as.data.frame(mcols(.ele))
            score <- apply(mc, 1, prod, na.rm = TRUE)
            mc <- mc * score
            ## change the signature if counts of 1 < counts of -1
            apply(mc, 2, function(.e){
                (sum(.e == -1)/length(.e))>0.5
                })
            })
        ## do sign switch
        for(chr in names(regions_sign_summary)){
            x <- regions_sign_summary[[chr]]
            for(samp in colnames(x)){
                if(x[samp]){
                    sample_data <- d[[samp]]
                    k <- seqnames(sample_data)==chr
                    mcols(sample_data[k])[, 'score'] <- -1 * mcols(sample_data[k])[, 'score']
                    d[[samp]] <- sample_data
                }
            }
        }
        ## export all A/B compartment scores
        null <- mapply(d, names(d), FUN=function(.d, .n){
            export(.d[!is.na(mcols(.d)[, 'score'])], file.path(OUTFOLDER, .n))
        })
    }
    """
}
