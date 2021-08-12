// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECKSUMS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    path "md5.*.txt", emit: md5

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    // Check md5
    def os_md5 = params.md5sum
    def os_md5_post = ""
    if (!params.md5sum){
        def osinfo = System.properties['os.name'].toLowerCase()
        os_md5 = osinfo.contains('mac')?'md5 -q':osinfo.contains('windows')?'certutil -hashfile':'md5sum'
        if(osinfo.contains('windows')){
            os_md5_post = 'MD5 | find /i /v `"md5`" | find /i /v `"certutil`"'
        }
    }
    """
    touch md5.${prefix}.txt
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    gunzip -c ${prefix}_1.fastq.gz > ${prefix}_1.fastq
    ${os_md5} ${prefix}_1.fastq ${os_md5_post} >>md5.${prefix}.txt
    gunzip -c ${prefix}_2.fastq.gz > ${prefix}_2.fastq
    ${os_md5} ${prefix}_2.fastq  ${os_md5_post} >>md5.${prefix}.txt
    if [ "${meta.md5_1}" != "null" ] && [ "${meta.md5_1}" != "" ]; then
        md5=(\$(${os_md5} ${prefix}_1.fastq.gz ${os_md5_post}))
        if [ "\$md5" != "${meta.md5_1}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_1}, but we got checksum \$md5!"
            exit 128
        fi
    fi
    if [ "${meta.md5_2}" != "null" ] && [ "${meta.md5_2}" != "" ]; then
        md5=(\$(${os_md5} ${prefix}_2.fastq.gz ${os_md5_post}))
        if [ "\$md5" != "${meta.md5_2}" ]
        then
            echo "${meta.id} has checksum ${meta.md5_2}, but we got checksum \$md5!"
            exit 128
        fi
    fi
    """
}
