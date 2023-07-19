process R1READS {
    tag "$meta.id"
    label 'process_low'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    gunzip -c $pair | \\
        awk 'BEGIN {OFS="\t"};  /^[^#]/ { print \$2, \$3, \$3+1, "*", "*", \$6}' | \\
        sort -k1,1 -k2,2n | \\
        gzip -nc > ${prefix}.R1.distal.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
    END_VERSIONS
    """
}
