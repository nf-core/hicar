process MERGE_READS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2' :
        'quay.io/biocontainers/samtools:1.10--h9402c20_2' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def sort_mem = task.memory * 0.8
    """
    gunzip -c ${bed} | \\
        sort -k1,1 -k2,2n --parallel $task.cpus -S ${sort_mem.toString().replaceAll(/ |B/, "")} | \\
        gzip -nc > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/[[:alpha:]|(|[:space:]]//g')
    END_VERSIONS
    """
}
