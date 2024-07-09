process HICEXPLORER_CHICVIEWPOINTBACKGROUNDMODEL {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(cool), path(anchors)

    output:
    tuple val(bin_size), path("${prefix}")      , emit: background_model
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}_background_model.txt"
    """
    chicViewpointBackgroundModel \\
        -t $task.cpus \\
        -m $cool \\
        -rp $anchors \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicViewpointBackgroundModel --version | sed 's/chicViewpointBackgroundModel //')
    END_VERSIONS
    """
}
