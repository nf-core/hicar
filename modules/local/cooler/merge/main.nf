process COOLER_MERGE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cooler=0.8.11"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*.cool"), emit: cool
    path "versions.yml"            , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}${meta.bin}"
    def args = task.ext.args ?: ''
    """
    cooler merge \\
        $args \\
        ${prefix}.cool \\
        ${cool}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(echo \$(cooler --version 2>&1) | sed 's/cooler, version //')
    END_VERSIONS
    """
}
