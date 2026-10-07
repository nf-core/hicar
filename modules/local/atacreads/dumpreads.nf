process DUMP_READS {
    tag "${meta.id}"
    label 'process_medium'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bed)
    val short_bed_postfix

    output:
    tuple val(meta), path("*.${short_bed_postfix}")   , emit: peak
    tuple val(meta), stdout                           , emit: counts
    path "versions.yml"                               , emit: versions

    script:
    def software = "awk"
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    gunzip -c $bed | \\
        awk -F "\t" 'BEGIN { OFS=FS } {print \$1,\$2,\$3 \\
        > "${prefix}."\$1".${short_bed_postfix}"}'
    echo \$(cat *.${short_bed_postfix} | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed -e "s/GNU Awk //g; s/, API.*\$//")
    END_VERSIONS
    """
}
