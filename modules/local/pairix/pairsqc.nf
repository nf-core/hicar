process PAIRSQC {
    tag "$meta.id"
    label 'process_low'
    errorStrategy { (task.exitStatus in 137..140 && task.attempt <= 3)  ? 'retry' : 'ignore' }

    conda "bioconda::pairix=0.3.7"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairix:0.3.7--py36h30a8e3e_3' :
        'biocontainers/pairix:0.3.7--py36h30a8e3e_3' }"

    input:
    tuple val(meta), path(pair), path(index)
    path chrom_sizes

    output:
    tuple val(meta), path("${meta.id}_report/*"), emit: qc
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    MAX_LOGDISTANCE=\$( cat ${chrom_sizes} | awk '{ sum += \$2 } END { printf "%.1f", log(sum)/log(10) }' )
    pairsqc.py \\
        -p $pair \\
        -c $chrom_sizes -t P \\
        -O $prefix \\
        -s $prefix \\
        -M \$MAX_LOGDISTANCE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairsqc: "0.2.2"
    END_VERSIONS
    """
}
