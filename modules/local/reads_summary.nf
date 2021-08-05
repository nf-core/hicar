// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READS_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.enable_conda ? "r::r-magrittr=1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-magrittr:1.5--r3.2.2_0"
    } else {
        container "quay.io/biocontainers/r-magrittr:1.5--r3.2.2_0"
    }

    input:
    path stat

    output:
    path "*.csv"                  , emit: summary
    path "*.version.txt"          , emit: version

    script:
    def software = "R"
    """
    read_summary.R

    echo \$(R --version 2>&1) | sed 's/R version //; s/Copyright.*\$//' > ${software}.version.txt
    """
}
