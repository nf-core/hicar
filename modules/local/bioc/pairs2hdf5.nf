process BIOC_PAIRS2HDF5 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bioconductor-trackviewer=1.28.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0' :
        'biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(pairs)
    path chromsizes

    output:
    tuple val(meta), path("*.h5") , emit: hdf5
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: 'keep-dup'
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Dec. 03, 2021 to convert pairs to hdf5
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    ## hdf5 format:
    ## - header, including total reads, chromosome name and sizes, tileWidth
    ##    * header/chrom_sizes COMPOUND
    ##    * header/header      STRING
    ##    * header/tile_width  INTEGER
    ##    * header/total       INTEGER
    ## - data, pairs in path data/chr1/chr2/tileIndex1_tileIndex2/
    ##    * position, in path data/chr1/chr2/tileIndex1_tileIndex2/position
    ##    * strand, in path data/chr1/chr2/tileIndex1_tileIndex2/strand
    #######################################################################
    #######################################################################

    pkgs <- c("rhdf5")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    comment_char <- "#"
    pattern <- "unselected.pairs.gz" ## this is from upstream output file.
    tileWidth <- 1e7 # this will create about 200 groups for human data
    block_size <- 1e7 # the size of block for tempfile
    keepDup <- FALSE # remove duplicates or not
    if(grepl("keep-dup", "$args")){
        keepDup <- TRUE
    }
    chrom_sizes <- read.delim("$chromsizes", header=FALSE)
    infs <- "$pairs"
    infs <- strsplit(infs, "\\\\s+")[[1]]

    getHeader <- function(inf){
        f <- gzfile(inf, open = "r")
        on.exit(close(f))
        header <- c()
        while(length(chunk <- readLines(f, n=1))){
            if(substr(chunk, 1, 1) == comment_char){
                header <- c(header, chunk)
            }else{
                break
            }
        }
        close(f)
        on.exit()
        header
    }
    getData <- function(f, block_size, n=7){
        if(n<0) return(data.frame())
        pc <- try({read.table(f, nrow=block_size, comment.char = "#",
                            colClasses=c(
                                "NULL", # reads name, skip
                                "character", # chrom1
                                "integer", # start1
                                "character", # chrom2
                                "integer", # start2
                                "character", #strand1
                                "character", #strand2
                                rep("NULL", n)))}, silent = TRUE)
        if(inherits(pc, "try-error")){
            data.frame()
        }else{
            pc
        }
    }
    getIndex <- function(pos, tileWidth){
        ceiling(pos/tileWidth)
    }
    getPath <- function(root, ...){
        paste(root, ..., sep="/")
    }
    createGroup <- function(obj, path){
        if(!H5Lexists(obj, path)){
            h5createGroup(obj, path)
        }
    }
    read_pair_write_tmp <- function(inf, out){
        ## check ncol
        h <- read.table(inf, nrow=1, comment.char = "#")
        n <- ncol(h) - 7
        f <- gzfile(inf, open = "r")
        on.exit({
            close(f)
        })
        filenames <- c()

        while(nrow(pc <- getData(f, block_size, n))>0){
            idx1 <- getIndex(pc[, 2], tileWidth)
            idx2 <- getIndex(pc[, 4], tileWidth)
            idx1_2 <- paste(pc[, 1], pc[, 3], paste(idx1, idx2, sep="_"), sep="___")
            pc <- split(pc[, -c(1, 3)], f=idx1_2)
            filenames <- unique(c(filenames, file.path(out, names(pc))))
            mapply(pc, names(pc), FUN=function(.data, .name){
                write.table(.data, file = file.path(out, .name),
                            append = TRUE, sep="\t", quote=FALSE,
                            col.names=FALSE, row.names=FALSE)
            })
        }

        close(f)
        on.exit()
        sort(filenames)
    }
    rewrite_hd5 <- function(filenames, out_h5, keepDup, root){
        total <- lapply(filenames, function(n){
            if(file.exists(n)){
                pc <- read.delim(n, header=FALSE)
                if(!keepDup) pc <- unique(pc)
                npc <- nrow(pc)
                chunk_size <- ceiling(sqrt(npc)/1000)*1000
                n <- getPath(root, gsub("___", "/", basename(n)))
                h5createGroup(out_h5, n)
                pos <- getPath(n, "position")
                if(npc>1000){
                    h5createDataset(out_h5, pos, dims = dim(pc[, c(1, 2)]),
                                    storage.mode = "integer",
                                    chunk = c(chunk_size, 2))
                }
                h5write(as.matrix(pc[, c(1, 2)]), out_h5, pos)
                strand <- getPath(n, "strand")
                if(npc>1000){
                    h5createDataset(out_h5, strand, dims = dim(pc[, c(3, 4)]),
                                    storage.mode = "character", size = 1,
                                    chunk = c(chunk_size, 2))
                }
                h5write(as.matrix(pc[, c(3, 4)]), out_h5, strand)
                npc
            }else{
                0
            }
        })
        sum(unlist(total))
    }
    root <- "data"
    for(inf in infs){
        out <- sub(pattern, "h5", basename(inf))
        if(!h5testFileLocking(dirname(inf))){
            h5disableFileLocking()
        }
        h5createFile(out)
        #header start as comment char'#'
        header <- getHeader(inf)
        h5createGroup(out, "header")
        h5write(header, out, "header/header")
        h5write(chrom_sizes, out, "header/chrom_sizes")
        h5write(as.integer(tileWidth), out, "header/tile_width")
        h5createGroup(out, root)
        # create groups for tiles
        n_chrom <- nrow(chrom_sizes)
        for(i in seq.int(n_chrom)){
            h5createGroup(out, getPath(root, chrom_sizes[i, 1]))
            for(j in seq.int(n_chrom)){
                h5createGroup(out, getPath(root, chrom_sizes[i, 1],
                                        chrom_sizes[j, 1]))
            }
        }

        # read pairs
        #columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
        tmp_dir <- "tmp_files"
        dir.create(tmp_dir)
        filenames <- try(read_pair_write_tmp(inf, tmp_dir))
        #rewrite
        total <- rewrite_hd5(filenames, out, keepDup, root)
        h5write(total, out, "header/total")
        h5closeAll()
        unlink(tmp_dir, recursive=TRUE)
    }
    """
}
