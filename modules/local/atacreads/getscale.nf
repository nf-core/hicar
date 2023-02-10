process GET_SCALE {
    tag "${meta.id}"
    label 'process_low'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), val(counts)
    val min_counts

    output:
    tuple val(meta), val(scale)  , emit: scale
    path "versions.yml"                               , emit: versions

    script:
    def software = "awk"
    scale = min_counts/counts.toInteger()
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    awk 'BEGIN { print $min_counts/$counts }'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed -e "s/GNU Awk //g; s/, API.*\$//")
    END_VERSIONS
    """
}
