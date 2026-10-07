process MAPS_MERGE {
    tag "$bin_size"
    label 'process_low'

    conda "pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'biocontainers/pandas:1.1.5' }"

    input:
    tuple val(bin_size), path(cut), path(mappability)

    output:
    tuple val(bin_size), path("${cut.getSimpleName()}")    , emit: map
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    merge_map.py \\
        -c $cut \\
        -m $mappability \\
        -o tmp.map
    awk "\\\$7>$args" tmp.map > ${cut.getSimpleName()}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
