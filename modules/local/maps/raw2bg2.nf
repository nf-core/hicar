process MAPS_RAW2BG2 {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::r-data.table=1.12.2"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-data.table:1.12.2' :
        'biocontainers/r-data.table:1.12.2' }"

    input:
    tuple val(meta), val(bin_size), path(peak)

    output:
    tuple val(meta), val(bin_size), path("${prefix}.bg2")          , emit: bg2
    tuple val(meta), val(bin_size), path("${prefix}.ginteractions"), emit: gi
    path "versions.yml"           , emit: versions

    script:
    def args   = task.ext.args ?: ''
    prefix   = task.ext.prefix ? "${task.ext.prefix}" : "${meta.id}_${bin_size}"
    """
    #!/usr/bin/env Rscript

    #########################################
    # Author: jianhong ou
    # create single files from reg_raw to bedgraph for cool load and juicerbox
    # the signal is the ratio2 (matrix count / matrix fit model expected2)
    #########################################
    versions <- c("${task.process}:", "    MAPS: 1.1.0")
    writeLines(versions, "versions.yml") # write versions.yml

    options("scipen"=999)
    library(data.table)
    RESOLUTION = as.numeric($bin_size)
    BG2_OUT = "${prefix}.bg2"
    GI_OUT  = "${prefix}.ginteractions"
    infs = strsplit("${peak}", "\\\\s+")[[1]]

    peaks_final_out <- lapply(infs, function(inf){
        peaks = as.data.table(read.table(inf, header=TRUE, stringsAsFactors=FALSE))
        if(nrow(peaks)==0){
            peaks\$bin1_end <- peaks\$bin2_end <- numeric(0)
        }else{
            peaks\$bin1_end = peaks\$bin1_mid + RESOLUTION
            peaks\$bin2_end = peaks\$bin2_mid + RESOLUTION
        }
        peaks_final = subset(peaks, select = c( "chr", "bin1_mid", "bin1_end",
                                                "chr", "bin2_mid", "bin2_end",
                                                "ratio2"))
        colnames(peaks_final) = c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count')
        peaks_final
    })
    peaks_final_out <- do.call(rbind, peaks_final_out)
    peaks_final_out <- peaks_final_out[peaks_final_out\$count>0, , drop=FALSE]
    peaks_final_out <- peaks_final_out[order(peaks_final_out\$chrom1,
                                            peaks_final_out\$start1,
                                            peaks_final_out\$chrom2,
                                            peaks_final_out\$start2), , drop=FALSE]
    write.table(peaks_final_out, BG2_OUT, row.names = FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    #print <strand1> <chr1> <pos1> <frag1> <strand2> <chr2> <pos2> <frag2> <score>
    peaks_final_out = cbind(0, peaks_final_out[, c('chrom1', 'start1')], 0,
                            0, peaks_final_out[, c('chrom2', 'start2')], 1,
                            peaks_final_out[, 'count'])
    write.table(peaks_final_out, GI_OUT, row.names = FALSE, col.names = FALSE, quote=FALSE, sep='\t')
    """
}
