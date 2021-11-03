// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCOS_PREPARE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bioconductor-rtracklayer=1.50.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    } else {
        container "quay.io/biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    }

    input:
    tuple val(meta), path(bedpe), val(ucscname), path(gtf), path(chromsize)

    output:
    tuple val(meta), path("circos/*")               , emit: circos
    path "*.version.txt"                            , emit: version

    script:
    """
    circos.r \\
        -i $bedpe \\
        -g $gtf \\
        -c $chromsize \\
        -u $ucscname
    """
}
