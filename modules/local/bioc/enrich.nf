// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_ENRICH {
    tag "$bin_size"
    label 'process_high'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:bin_size) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bioconductor-clusterprofiler=3.18.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    }

    input:
    tuple val(bin_size), path(diff)

    output:
    tuple val(bin_size), path("diffhic_bin${bin_size}/enrichment/*"), emit: enrichment
    path "*.version.txt"                                            , emit: version

    script:
    """
    install_packages.r clusterProfiler pathview biomaRt optparse
    enrich.r -s ${params.species} -o "diffhic_bin${bin_size}/enrichment" $options.args

    # *.version.txt
    """
}
