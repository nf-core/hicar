// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECKSUMS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "python=3.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8"
    } else {
        container "quay.io/biocontainers/python:3.8"
    }

    input:
    tuple val(meta), path(reads)

    output:
    path "md5.*.txt"              , emit: md5
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    touch md5.${prefix}.txt
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    gunzip -c ${prefix}_1.fastq.gz > ${prefix}_1.fastq
    md5.py ${prefix}_1.fastq >>md5.${prefix}.txt
    gunzip -c ${prefix}_2.fastq.gz > ${prefix}_2.fastq
    md5.py ${prefix}_2.fastq >>md5.${prefix}.txt
    if [ "${meta.md5_1}" != "null" ] && [ "${meta.md5_1}" != "" ]; then
        md5=(\$(md5.py ${prefix}_1.fastq.gz))
        if [ "\$md5" != "${meta.md5_1}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_1}, but we got checksum \$md5!"
            exit 128
        fi
    fi
    if [ "${meta.md5_2}" != "null" ] && [ "${meta.md5_2}" != "" ]; then
        md5=(\$(md5.py ${prefix}_2.fastq.gz))
        if [ "\$md5" != "${meta.md5_2}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_2}, but we got checksum \$md5!"
            exit 128
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS
    """
}
