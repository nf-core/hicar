// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READS_SUMMARY {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "r::r-magrittr=1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-magrittr:1.5--r3.2.2_0"
    } else {
        container "quay.io/biocontainers/r-magrittr:1.5--r3.2.2_0"
    }

    input:
    path stat

    output:
    path "*.{csv,json}"           , emit: summary
    path "versions.yml"           , emit: versions

    script:
    """
    read_summary.R

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        R: \$(echo \$(R --version 2>&1) | sed 's/R version //; s/Copyright.*\$//')
    END_VERSIONS
    """
}
