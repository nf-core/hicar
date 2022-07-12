process COOLER_BALANCE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: cool
    path "versions.yml"                            , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}${meta.bin}_balanced"
    def args = task.ext.args ?: ''
    extension = cool.getExtension()
    """
    cp ${cool} ${prefix}.${extension}
    cooler balance \\
        --nproc $task.cpus \\
        $args \\
        ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(echo \$(cooler --version 2>&1) | sed 's/cooler, version //')
    END_VERSIONS
    """
}
