process HICEXPLORER_HICCORRECTMATRIX {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.h5")               , emit:corrected
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}_corrected"
    """
    hicCorrectMatrix correct \\
        -m $cool \\
        $args \\
        -o ${prefix}.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicCorrectMatrix --version 2>&1 | sed 's/hicCorrectMatrix //')
    END_VERSIONS
    """
}
