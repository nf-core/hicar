// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process ENSEMBL_UCSC_CONVERT {
    tag "$fname"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bioconductor-rtracklayer=1.50.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    } else {
        container "quay.io/biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    }

    input:
    tuple val(bin_size), path(fname)

    output:
    tuple val(bin_size), path("{UCSC,ensembl}.${fname}"), emit: tab
    path "*.version.txt"          , emit: version

    script:
    """
    seqlevels_convert.r \\
        $options.args \\
        $fname
    """
}
