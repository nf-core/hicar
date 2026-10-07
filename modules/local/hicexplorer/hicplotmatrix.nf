process HICEXPLORER_HICPLOTMATRIX {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool), path(additional_files)

    output:
    tuple val(meta), path("*.{png,jpg,jpeg}")   , emit:plots
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    hicPlotMatrix \\
        -m $cool \\
        $args \\
        -o ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotMatrix --version 2>&1 | sed 's/hicPlotMatrix //')
    END_VERSIONS
    """
}
