// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process DIFFHICAR {
    tag "$bin_size"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bioconductor-edger=3.32.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-edger:3.32.1--r40h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-edger:3.32.1--r40h399db7b_0"
    }

    input:
    tuple val(bin_size), path(peaks, stageAs: "peaks/*"), path(long_bedpe, stageAs: "long/*")

    output:
    tuple val(bin_size), path("diffhic_bin${bin_size}/*"), emit: diff
    path "*.version.txt"                  , emit: version

    script:
    """
    install_packages.r edgeR
    diffhicar.r diffhic_bin${bin_size} \\
        $options.args
    ## must output the packages version as *.version.txt
    """
}
