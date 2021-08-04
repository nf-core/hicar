// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRSQC {
    tag "$meta.id"
    label 'process_low'
    errorStrategy { (task.exitStatus in 137..140 && task.attempt <= 3)  ? 'retry' : 'ignore' }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::pairix=0.3.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pairix:0.3.7--py36h30a8e3e_3"
    } else {
        container "quay.io/biocontainers/pairix:0.3.7--py36h30a8e3e_3"
    }

    input:
    tuple val(meta), path(pair), path(index)
    path chrom_sizes

    output:
    tuple val(meta), path("${meta.id}_report/*"), emit: qc
    path "*.version.txt"          , emit: version

    script:
    def software = "pairsqc"
    """
    MAX_LOGDISTANCE=`cat genome.fa.sizes | awk '{ sum += \$2 } END { printf "%.1f", log(sum)/log(10) }'`
    pairsqc.py \\
        -p $pair \\
        -c $chrom_sizes -t P \\
        -O ${meta.id} \\
        -s ${meta.id} \\
        -M \$MAX_LOGDISTANCE

    echo "0.2.2" > ${software}.version.txt
    """
}
