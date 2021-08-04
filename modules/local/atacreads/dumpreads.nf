// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process DUMPREADS {
    tag "${meta.id}"
    label 'process_medium'
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
    tuple val(meta), path("${meta.id}-${params.species}/*.shrt.vip.bed") , emit: peak
    path  "*.version.txt"                                                , emit: version

    script:
    def software = "awk"
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def outdir   = "${meta.id}-${params.species}"
    """
    mkdir -p $outdir
    cat $bed | zcat | awk -v setname=${prefix} -v outdir=${outdir}  -F \$"\t"  'BEGIN {OFS=FS};{print \$1,\$2,\$3 > outdir"/"setname"."\$1".shrt.vip.bed"}'

    echo \$(awk --version 2>&1) | sed -e "s/GNU Awk //g; s/, API.*\$//" > ${software}.version.txt
    """
}
