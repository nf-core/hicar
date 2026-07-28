process BIOC_ATACSEQTFEA {
    label 'process_high_memory'
    label 'process_long'
    label 'process_single'
    label 'error_ignore'

    conda "bioconda::bioconductor-atacseqtfea=1.0.1"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-atacseqtfea:1.0.1--r42hdfd78af_0' :
        'biocontainers/bioconductor-atacseqtfea:1.0.1--r42hdfd78af_0' }"

    input:
    tuple val(meta), path(peak), path(bam)
    tuple val(genome), path(fasta), val(gtf)

    output:
    tuple val(meta), path("${prefix}/")                         , emit: tfea
    path  "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Aug. 24, 2021 for ATACseqTFEA
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################

    pkgs <- c("ATACseqTFEA", "rtracklayer", "BSgenome", "Rsamtools")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    # Options
    PVALUE <- 0.05
    motifpath <- "best_curated_Human.rds"
    gtf <- "${gtf}"
    output <- "${prefix}"
    genome <- "$genome"
    PEAKS <- strsplit("$peak", " ")[[1]]
    bams <- dir(pattern="bam\$")
    stopifnot(length(bams)>=1)
    args <- strsplit("${args}", "\\\\s+")[[1]]
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
    option_list <- list("motifpath"=c("--motifpath", "-m", "character"),
                        "seqlevels"=c("--seqlevels", "-l", "character"))
    opt <- parse_args(option_list, args)

    if(!is.null(opt\$motifpath)){
        motifpath <- opt\$motifpath
    }
    if(!is.null(opt\$seqlevels)){
        seqlevels <- eval(opt\$seqlevels)
    }

    if(file.exists(system.file("extdata", motifpath, package = "ATACseqTFEA"))){
        motifpath <- system.file("extdata", motifpath, package = "ATACseqTFEA")
    }else{
        motifpath <- motifpath
    }

    ## prepare binding sites
    motifs <- readRDS(motifpath)
    gtf <- import(gtf)
    peaks <- lapply(PEAKS, import)
    peaks <- reduce(unlist(GRangesList(peaks)))
    genome <- tryCatch({
        getBSgenome(genome)
    }, error=function(e){
        message("Try to forge a genome")
        message("It will fail if the BSgenome package is not proper installed, check this thread:https://support.bioconductor.org/p/124169/")
        tmp <- paste0(genome, '.seeds')
        pkgName <- paste0("Package: BSgenome.", genome)
        writeLines(c(
            pkgName,
            "Title: Full genome sequences for temp usage",
            "Description: Full genome sequences for TFEA",
            "Version: 0.0.1",
            paste("BSgenomeObjname:", genome),
            "SrcDataFiles: $fasta",
            "seqs_srcdir: .",
            paste("genome:", genome),
            paste("organism:", genome),
            paste("common_name:", genome),
            paste("organism_biocview:", genome),
            "provider: UNKONWN",
            paste("release_date:", Sys.Date()),
            paste0("seqnames: c('", paste(seqlevels(gtf), collapse = "','"), "')")
        ), tmp)
        # try to touch the inst/extdata folder in the template
        tryCatch({
            bsgenome_path <- system.file(package="BSgenome")
            dir.create(file.path(bsgenome_path, 'pkgtemplates', 'BSgenome_datapkg', 'inst', 'extdata'), recursive=TRUE)
        }, error = function(.e){
            message(.e)
        })
        forgeBSgenomeDataPkg(tmp)
        pkgName <- paste0("BSgenome.", genome)
        libpath <- 'localLib'
        dir.create(libpath)
        install.packages(pkgName, repos = NULL, type="source", lib=libpath)
        library(pkgName, character.only = TRUE, lib.loc=libpath)
        get(pkgName)
    })
    bindingSites <-
        prepareBindingSites(
            motifs,
            genome,
            grange=peaks,
            p.cutoff = PVALUE
        )
    ## prepare R2 read bam file
    bams <- sapply(bams, function(.ele) sortBam(.ele, sub(".bam\$", ".srt", .ele)))
    indexBam(bams)
    bams <- lapply(bams, function(.ele) {
        filterBam(.ele, sub(".bam\$", ".R2.bam", .ele),
            param=ScanBamParam(
                flag = scanBamFlag(isSecondMateRead=TRUE),
                what=scanBamWhat()))
    })
    bams <- unlist(bams)
    group <- sub("_REP.*\$", "", bams)
    stopifnot(length(unique(group))>=1)
    group <- factor(group)
    ## TFEA
    if(levels(group)>=2){
        contrasts <- combn(levels(group), m=2, simplify=FALSE)
    }else{
        contrasts <- levels(group)
    }
    for(cont in contrasts){
        gp1 <- cont[1]
        bamExp <- bams[group %in% gp1]
        if(length(cont)==2){
            gp2 <- cont[2]
            bamCtl <- bams[group %in% gp2]
        }else{
            gp2 <- "NA"
            bamCtl <- NULL
        }
        res <- TFEA(
            bamExp=bamExp,
            bamCtl=bamCtl,
            bindingSites = bindingSites
        )
        ## save results
        output_cont <- file.path(output, paste0(gp1, ".vs.", gp2))
        dir.create(output_cont, recursive = TRUE)
        saveRDS(res, file.path(output_cont, "ATACseqTFEA_results.rds"))
        write.csv(
            res[["resultsTable"]],
            file.path(output_cont, "ATACseqTFEA_results.csv"),
            row.names = FALSE)
        TFS <- res\$resultsTable\$TF[res\$resultsTable\$adjPval<PVALUE]
        pdf(file.path(output_cont, "ATACseqTFEA_ESvolcanoplot.pdf"))
        ESvolcanoplot(TFEAresults=res, TFnameToShow=TFS)
        dev.off()
        outfolder <- file.path(output_cont, "ESplots")
        for(TF in TFS){
            plotES(res, TF=TF, outfolder=outfolder)
        }
    }
    """
}
