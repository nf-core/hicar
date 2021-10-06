// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_CHIPPEAKANNO {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) }

    conda (params.enable_conda ? "bioconda::bioconductor-chippeakanno=3.26.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-chippeakanno:3.26.0--r41hdfd78af_0"
    }

    input:
    tuple val(bin_size), path(diff)
    path gtf

    output:
    tuple val(bin_size), path("${prefix}/anno/*"), emit: anno
    tuple val(bin_size), path("${prefix}/anno/**.anno.csv"), emit: csv
    path "${prefix}/anno/*.png", optional:true, emit: png
    path "*.version.txt"               , emit: version

    script:
    prefix   = options.suffix ? "${options.suffix}${bin_size}" : "diffhic_bin${bin_size}"
    """
    install_packages.r ChIPpeakAnno ggplot2
    annopeaks.r ${gtf} ${prefix}
    """
}
