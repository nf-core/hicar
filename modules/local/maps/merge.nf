process MAPS_MERGE {
    tag "$bin_size"
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(bin_size), path(cut), path(mappability)
    path merge_map_py_source

    output:
    tuple val(bin_size), path("${cut.getSimpleName()}")    , emit: map
    path "versions.yml"                                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    python $merge_map_py_source \\
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
