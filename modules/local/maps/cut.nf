process MAPS_CUT {
    tag "$bin_size"
    label 'process_high_cpus'
    label 'process_long'
    label 'error_retry'

    conda "conda-forge::biopython=1.70"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_1' :
        'biocontainers/biopython:1.70--np112py36_1' }"

    input:
    tuple val(bin_size), val(site), path(fasta)
    val enzyme

    output:
    tuple val(bin_size), path('*.cut')        , emit: cut
    path "versions.yml"                       , emit: versions

    script:
    """
    cut=($site)
    restriction_cut_multipleenzyme.py \\
        -f ${fasta} \\
        -s \${cut[0]} \\
        -p \${cut[1]} \\
        -b ${bin_size} \\
        -o ${bin_size}_${enzyme}.cut \\
        -c $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
