process CHECKSUMS {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::coreutils=8.31" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0' :
        'quay.io/biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "*.md5"                  , emit: md5
    path "versions.yml"           , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    gunzip -c ${prefix}_1.fastq.gz > ${prefix}_1.fastq
    gunzip -c ${prefix}_2.fastq.gz > ${prefix}_2.fastq
    md5sum ${prefix}_1.fastq ${prefix}_2.fastq > ${prefix}.md5
    rm ${prefix}_1.fastq
    rm ${prefix}_2.fastq
    if [ ! -z "${meta.md5_1 ?: ''}" ] && [ ! -z "${meta.md5_2 ?: ''}" ]; then
        cat <<-END_CHECKSUM | md5sum -c
    ${meta.md5_1}  ${prefix}_1.fastq.gz
    ${meta.md5_2}  ${prefix}_2.fastq.gz
    END_CHECKSUM
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS
    """
}
