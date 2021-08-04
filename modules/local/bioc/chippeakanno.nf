// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_CHIPPEAKANNO {
    tag "$bin_size"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bioconductor-chippeakanno=3.24.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.24.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-chippeakanno:3.24.1--r40hdfd78af_0"
    }

    input:
    tuple val(bin_size), path(diff)
    path gtf

    output:
    tuple val(bin_size), path("diffhic_bin${bin_size}/anno/*"), emit: anno
    path "*.version.txt"               , emit: version

    script:
    """
    install_packages.r ChIPpeakAnno ggplot2
    annopeaks.r ${gtf} diffhic_bin${bin_size}
    """
}
