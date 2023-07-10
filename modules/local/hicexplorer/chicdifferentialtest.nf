process HICEXPLORER_CHICDIFFERENTIALTEST {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(aggregate)

    output:
    tuple val(bin_size), path("${prefix}_differential.hdf5")    , emit: differential
    path "versions.yml"                                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}"
    """
    chicDifferentialTest \\
        $args \\
        --aggregatedFile $aggregate \\
        --outFileName ${prefix}_differential.hdf5 \\
        -t $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotViewpoint --version 2>&1 | sed 's/hicPlotViewpoint //')
    END_VERSIONS
    """
}
