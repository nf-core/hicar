process DUMP4HOMER {
    tag "${meta.id}"
    label 'process_medium'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path("*.HiCsummary.txt.gz")           , emit: hicsummary
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c $pair | \\
        awk -v F='\\t' -v OFS='\\t' \\
        '!/^[[:space:]]*#/ {print \$1,\$2,\$3,\$6,\$4,\$5,\$7}' | \\
        gzip -nc > ${prefix}.HiCsummary.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed -e "s/GNU Awk //g; s/, API.*\$//")
    END_VERSIONS
    """
}
