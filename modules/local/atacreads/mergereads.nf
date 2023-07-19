process MERGE_READS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.10"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2' :
        'biocontainers/samtools:1.10--h9402c20_2' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c ${bed} | \\
        sort -k1,1 -k2,2n | \\
        gzip -nc > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/[[:alpha:]|(|[:space:]]//g')
    END_VERSIONS
    """
}
