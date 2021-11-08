// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READS_STAT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "r::r-magrittr=1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-magrittr:1.5--r3.2.2_0"
    } else {
        container "quay.io/biocontainers/r-magrittr:1.5--r3.2.2_0"
    }

    input:
    tuple val(meta), path(raw), path(dedup)

    output:
    tuple val(meta), path("*.csv"), emit: stat
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    read_stat.R $raw $dedup ${prefix}.reads_stats.csv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        R: \$(echo \$(R --version 2>&1) | sed 's/R version //; s/Copyright.*\$//')
    END_VERSIONS
    """
}
