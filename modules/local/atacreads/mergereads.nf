// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MERGEREADS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "*.version.txt"             , emit: version

    script:
    def software = "gzip"
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    cat ${bed} | zcat | \\
    sort -k1,1 -k2,2n |  gzip -nc > ${prefix}.merged.ATAC.bed.gz

    echo \$(gzip --version 2>&1) | sed 's/[[:alpha:]|(|[:space:]]//g' > ${software}.version.txt
    """
}
