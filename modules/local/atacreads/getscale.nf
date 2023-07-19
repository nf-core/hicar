process GET_SCALE {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), val(counts)
    val min_counts

    output:
    tuple val(meta), val(scale)  , emit: scale
    path "versions.yml"          , emit: versions

    script:
    scale = min_counts/counts.toInteger()
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextflow: "${nextflow.version}"
    END_VERSIONS
    """
}
