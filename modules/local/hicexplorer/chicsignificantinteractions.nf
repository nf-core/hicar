process HICEXPLORER_CHICSIGNIFICANTINTERACTIONS {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(background), path(interactions)

    output:
    path "${prefix}_significant.hdf5"            , emit: significant
    path "${prefix}_target.hdf5"                 , emit: target
    tuple val(bin_size), path(interactions), path("${prefix}_target.hdf5"), emit: interaction_target
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${bin_size}"
    """
    chicSignificantInteractions \\
        --interactionFile $interactions \\
        -bmf $background \\
        --outFileNameSignificant ${prefix}_significant.hdf5 \\
        --outFileNameTarget ${prefix}_target.hdf5 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicSignificantInteractions --version 2>&1 | sed 's/chicSignificantInteractions //')
    END_VERSIONS
    """
}
