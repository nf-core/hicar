process HICEXPLORER_CHICQUALITYCONTROL {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(cool), path(anchors)

    output:
    tuple val(bin_size), path(cool), path("${prefix}"), emit: referencepoints
    path "versions.yml"                               , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}_referencepoints.bed"
    """
    chicQualityControl \\
        -t $task.cpus \\
        -m $cool \\
        -rp $anchors \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicQualityControl --version 2>&1 | sed 's/chicQualityControl //')
    END_VERSIONS
    """
}
