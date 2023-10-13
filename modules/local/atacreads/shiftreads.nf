process SHIFT_READS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::coreutils=8.31"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0' :
        'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(pair)
    val do_shift

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val(meta), stdout          , emit: counts
    path "versions.yml"              , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    def buffer   = task.memory.toGiga().intdiv(2)
    def command  = do_shift ?
        'BEGIN {OFS="\t"};  /^[^#]/ { if (\$7 == "+") {\$5 = \$5 + 4} else if (\$7 == "-") {\$5 = \$5 - 5};  print \$4, \$5, \$5+1, "*", "0", \$7}' :
        'BEGIN {OFS="\t"};  /^[^#]/ { print \$4, \$5, \$5+1, "*", "0", \$7 }'
    """
    gunzip -c $pair | \\
        awk '$command' | \\
        sort --parallel=$task.cpus --buffer-size=${buffer}G -k1,1 -k2,2n | \\
        uniq | \\
        gzip -nc > ${prefix}.R2.ATAC.bed.gz
    echo \$(cat $pair | wc -l)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
    END_VERSIONS
    """
}
