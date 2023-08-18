process HICEXPLORER_CHICAGGREGATESTATISTIC {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(interactions), path(target)

    output:
    tuple val(bin_size), path("${prefix}_aggregate.hdf5")    , emit: aggregate
    path "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}"
    """
    chicAggregateStatistic \\
        --interactionFile $interactions \\
        --targetFile $target \\
        --outFileName ${prefix}_aggregate.hdf5 \\
        -t $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicAggregateStatistic --version 2>&1 | sed 's/chicAggregateStatistic //')
    END_VERSIONS
    """
}
