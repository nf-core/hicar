process SHIFT_READS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(pair)
    val do_shift

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    def command  = do_shift ?
        'BEGIN {OFS="\t"};  /^[^#]/ { if (\$7 == "+") {\$5 = \$5 + 4} else if (\$7 == "-") {\$5 = \$5 - 5};  print \$4, \$5, \$5+1, "*", "0", \$7}' :
        'BEGIN {OFS="\t"};  /^[^#]/ { print \$4, \$5, \$5+1, "*", "0", \$7 }'
    def sort_mem = task.memory * 0.8
    """
    gunzip -c $pair | \\
        awk '$command' | \\
        sort -k1,1 -k2,2n -S ${sort_mem.toString().replaceAll(/ |B/, "")} | \\
        uniq | \\
        gzip -nc > ${prefix}.R2.ATAC.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
    END_VERSIONS
    """
}
