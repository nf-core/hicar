process HICEXPLORER_HICTRANSFORM {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*${prefix}.h5")                    , emit:transformed
    tuple val(meta), path("obs_exp_${prefix}.h5")             , emit:obs_exp, optional: true
    tuple val(meta), path("obs_exp_lieberman_${prefix}.h5")   , emit:obs_exp_lieberman, optional: true
    tuple val(meta), path("obs_exp_non_zero_${prefix}.h5")    , emit:obs_exp_non_zero, optional: true
    tuple val(meta), path("pearson_${prefix}.h5")             , emit:pearson, optional: true
    tuple val(meta), path("covariance_${prefix}.h5")          , emit:covariance, optional: true
    path("versions.yml")                                      , emit:versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    ## step1 create bigwig files
    hicTransform \\
        -m $cool \\
        $args \\
        --outFileName ${prefix}.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicTransform --version 2>&1 | sed 's/hicTransform //')
    END_VERSIONS
    """
}
