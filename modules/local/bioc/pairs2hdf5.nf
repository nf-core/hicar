process BIOC_PAIRS2HDF5 {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(pairs)
    path chromsizes

    output:
    tuple val(meta), path("*.h5") , emit: hdf5
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: 'keep-dup'
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on Dec. 03, 2021 to convert pairs to hdf5
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
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
    block_size <- 1e6 # the size of block for tempfile
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
    getDataIndex <- function(block_idx, data, subG){
        n <- nrow(data)-1
        #start index for next data
        idx <- ifelse(is.null(block_idx[[subG]]),
                        0, block_idx[[subG]][["idx"]])
        region <- c(start=idx, end=idx+n)
        # block_index
        block_index <- floor(region/block_size)
        if(block_index[1]==block_index[2]){
            block <- list(region-block_index*block_size+1)
            names(block) <- block_index[1]
        }else{
            block <- list(c(region[1]-block_index[1]*block_size+1,
                            block_size),
                        c(1, region[2]-block_index[2]*block_size+1))
            names(block) <- block_index
        }
        list(idx=region[2]+1, block=block, block_idx=block_index[2])
    }
    createGroup <- function(obj, path){
        if(!H5Lexists(obj, path)){
            h5createGroup(obj, path)
        }
    }
    read_pair_write_hd5 <- function(inf, out, root){
        ## check ncol
        h <- read.table(inf, nrow=1, comment.char = "#")
        n <- ncol(h) - 7
        f <- gzfile(inf, open = "r")
        out <- H5Fopen(out)
        on.exit({
            close(f)
            H5Fclose(out)
        })
        block_idx <- list()
        createGroup(out, root)

        while(nrow(pc <- getData(f, block_size, n))>0){
            pc <- split(pc, pc[, 1])
            for(i in names(pc)){
                createGroup(out, getPath(root, i))
                pc_j <- split(pc[[i]], pc[[i]][, 3])
                for(j in names(pc_j)){
                    createGroup(out, getPath(root, i, j))
                    idx1 <- getIndex(pc_j[[j]][, 2], tileWidth)
                    idx2 <- getIndex(pc_j[[j]][, 4], tileWidth)
                    pc_sub <- split(pc_j[[j]][, -c(1, 3)],
                                    paste(idx1, idx2, sep="_"))
                    for(k in names(pc_sub)){
                        subG <- getPath(root, i, j, k)
                        createGroup(out, subG)
                        idx <- getDataIndex(block_idx, pc_sub[[k]], subG)
                        start <- as.matrix(pc_sub[[k]][, c(1, 2)])
                        strand <- as.matrix(pc_sub[[k]][, c(3, 4)])
                        block <- idx[["block"]]
                        for(b in names(block)){
                            subG_block <- getPath(subG, b)
                            createGroup(out, subG_block)
                            subG_start <- getPath(subG_block, "start")
                            subG_strand <- getPath(subG_block, "strand")
                            if(block[[b]][1]==1){
                                h5createDataset(out, subG_start, dims=c(block[[b]][2], 2),
                                                maxdims=c(block_size, 2),
                                                storage.mode="integer")
                                h5createDataset(out, subG_strand, dims=c(block[[b]][2], 2),
                                                maxdims=c(block_size, 2),
                                                storage.mode="character", size=1)
                            }else{
                                h5set_extent(out, subG_start, dims=c(block[[b]][2], 2))
                                h5set_extent(out, subG_strand, dims=c(block[[b]][2], 2))
                            }
                            keep <- seq.int(diff(block[[b]])+1)
                            h5write(start[keep, , drop=FALSE], out, subG_start,
                                    index=list(seq(block[[b]][1], block[[b]][2]), NULL))
                            h5write(strand[keep, , drop=FALSE], out, subG_strand,
                                    index=list(seq(block[[b]][1], block[[b]][2]), NULL))
                            start <- start[-keep, , drop=FALSE]
                            strand <- strand[-keep, , drop=FALSE]
                        }
                        block_idx[[subG]] <-  idx
                    }
                }
            }

        }
        close(f)
        H5Fclose(out)
        on.exit()
        block_idx
    }
    read_hd5_write_hd5 <- function(in_h5, out_h5, block_idx, keepDup){
        input <- H5Fopen(in_h5, flags="H5F_ACC_RDONLY")
        on.exit(H5Fclose(input))
        total <- lapply(names(block_idx), function(n){
            if(H5Lexists(input, n)){
                pc <- lapply(seq.int(block_idx[[n]][["block_idx"]]+1)-1, function(block){
                    pos <- h5read(input, getPath(n, block, "start"))
                    strand <- h5read(input, getPath(n, block, "strand"))
                    cbind(as.data.frame(pos), as.data.frame(strand))
                })
                pc <- do.call(rbind, pc)
                if(!keepDup) pc <- unique(pc)
                npc <- nrow(pc)
                chunk_size <- ceiling(sqrt(npc)/1000)*1000
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
        tmp_h5 <- tempfile(tmpdir=getwd(), fileext = ".h5")
        h5createFile(tmp_h5)
        block_idx <- try(read_pair_write_hd5(inf, tmp_h5, root))
        #rewrite
        total <- read_hd5_write_hd5(tmp_h5, out, block_idx, keepDup)
        h5write(total, out, "header/total")
        h5closeAll()
        unlink(tmp_h5)
    }
    """
}
