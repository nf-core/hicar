process HICEXPLORER_CHICVIEWPOINT {
    tag "${bin_size}"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(cool), path(anchors), path(background)

    output:
    tuple val(bin_size), path(background), path("${prefix}")  , emit: interactions
    path "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}_interactions.hdf5"
    """
    chicViewpoint \\
        -t $task.cpus \\
        -m $cool \\
        -rp $anchors \\
        -bmf $background \\
        --outFileName $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotViewpoint --version 2>&1 | sed 's/hicPlotViewpoint //')
    END_VERSIONS
    """
}
