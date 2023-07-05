process HICEXPLORER_CHICPLOTVIEWPOINT {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(background), path(interactions), path(differential)

    output:
    path "${prefix}"              , emit: v4c
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}_v4c.tar.gz"
    """
    chicPlotViewpoint \\
        $args \\
        --interactionFile $interactions \\
        --differentialTestResult $differential \\
        --backgroundModelFile $background \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicPlotViewpoint --version 2>&1 | sed 's/chicPlotViewpoint //')
    END_VERSIONS
    """
}
