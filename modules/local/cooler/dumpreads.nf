process DUMPREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bedpe)
    val long_bedpe_postfix

    output:
    tuple val(meta), path("*.${long_bedpe_postfix}")       , emit: bedpe
    tuple val(meta), path("*.ginteractions")               , emit: gi
    path  "versions.yml"                                   , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    awk -F "\\t" \\
        '! /^chrom1/ {print > "${prefix}."\$1"_"\$4".${long_bedpe_postfix}"}' \\
        $bedpe

    awk -F "\\t" '{print 0, \$1, \$2, 0, 0, \$4, \$5, 1, \$7}' $bedpe > \\
        ${prefix}.${meta.bin}.ginteractions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed -e "s/GNU Awk //g; s/, API.*\$//")
    END_VERSIONS
    """
}
