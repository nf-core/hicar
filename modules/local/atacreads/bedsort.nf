process BEDFILES_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0' :
        'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(intervals)
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    path  "versions.yml"                   , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def buffer   = task.memory.toGiga().intdiv(2)
    """
    ## ref: https://www.biostars.org/p/66927/
    LC_ALL=C sort \\
        --parallel=$task.cpus \\
        --buffer-size=${buffer}G \\
        -k1,1 -k2,2n \\
        $intervals \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(sort --version | tr '\\n' ' ' | sed -e "s/^[^0-9]*//; s/ Copyright.*\$//")
    END_VERSIONS
    """
}
