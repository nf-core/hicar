process MAPS_REFORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::r-data.table=1.12.2"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-data.table:1.12.2' :
        'biocontainers/r-data.table:1.12.2' }"

    input:
    tuple val(meta), val(bin_size), path(peak)
    val maps_3d_ext

    output:
    tuple val(meta), val(bin_size), path("*.bedpe")             , emit: bedpe
    tuple val(meta), val(bin_size), path("*.${maps_3d_ext}")    , emit: peaks
    path "versions.yml"                                         , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    #########################################
    # Author: [Ivan Juric](https://github.com/ijuric)
    # File: MAPS_peak_formatting.r
    # Source: https://github.com/ijuric/MAPS/blob/master/bin/MAPS/MAPS_peak_formatting.r
    # Source+commit: https://raw.githubusercontent.com/ijuric/MAPS/46f92a6c7965a3b855ba7558c50d6c021b1d677c/bin/MAPS/MAPS_peak_formatting.r
    # Data: 11/08/2021, commit:46f92a6
    # modified by Jianhong:
    # ## 1. set the Rscript environment.
    # ## 2. simplify the input and output parameters by input filename and output filename
    # ## 3. handle multiple input files
    # ## 4. handle the empty input file
    #########################################
    versions <- c("${task.process}:", "    MAPS: 1.1.0")
    writeLines(versions, "versions.yml") # write versions.yml

    options("scipen"=999)
    library(data.table)
    RESOLUTION = as.numeric($bin_size)
    infs = strsplit("${peak}", "\\\\s+")[[1]]

    for(inf in infs){
        peaks_raw = read.table(inf,header=TRUE, stringsAsFactors=FALSE)

        peaks = as.data.table(subset(peaks_raw, ClusterType != 'Singleton' | (ClusterType == 'Singleton' & fdr < 1e-4))) # remove singletons
        if(nrow(peaks)==0){
            peaks\$bin1_end <- peaks\$bin2_end <- peaks\$summit <- numeric(0)
        }else{
            peaks[, summit := 1*(fdr == min(fdr)), by = lab]
            peaks\$summit[ peaks\$ClusterType == 'Singleton'] = 1

            singleton_labs = peaks\$lab[ peaks\$ClusterType == 'Singleton']
            peaks\$lab[ peaks\$ClusterType == 'Singleton'] = paste(singleton_labs, 1:length(singleton_labs),sep='')

            peaks\$bin1_end = peaks\$bin1_mid + RESOLUTION
            peaks\$bin2_end = peaks\$bin2_mid + RESOLUTION
        }
        peaks_final = subset(peaks, select = c("chr", "bin1_mid", "bin1_end", "chr", "bin2_mid", "bin2_end", "count", "expected2", "fdr", "lab", "ClusterSize", "ClusterType", "NegLog10P", "summit"))
        colnames(peaks_final) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'count', 'expected', 'fdr', 'ClusterLabel', 'ClusterSize', 'ClusterType', 'ClusterNegLog10P', 'ClusterSummit')
        fout = sub(".peaks",paste0('.', "$maps_3d_ext"),inf)
        write.table(peaks_final, fout, row.names = FALSE, col.names = TRUE, quote=FALSE, sep='\t')
        peaks_final\$name='.'
        peaks_final <- peaks_final[, c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'ClusterNegLog10P')]
        write.table(peaks_final,
            sub(".peaks", ".bedpe", inf),
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    }
    """
}
