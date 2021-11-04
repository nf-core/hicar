// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIR2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(peak), path(pairs)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai")    , emit: bam
    path "*.version.txt"                                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pair2bam.r $peak

    # *.version.txt files will be created in the rscripts
    """
}
