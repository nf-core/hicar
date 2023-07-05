process RECENTER_PEAK {
    label 'process_low'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    path peak

    output:
    path "recentered_peak.bed"    , emit: peak
    path "versions.yml"           , emit: versions

    script:
    """
    awk 'BEGIN {OFS="\\t"} {print \$1,int((\$2+\$3)/2),int((\$2+\$3)/2)+1,\$4}' $peak > recentered_peak.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
    END_VERSIONS
    """
}
