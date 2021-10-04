// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALL_HIPEAK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-monocle=2.20.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-monocle:2.20.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-monocle:2.20.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(r2peak, stageAs:"R2peak/*"), path(r1peak, stageAs: "R1peak/*"), path(distalpair, stageAs: "pairs/*"), path(mappability), path(fasta), path(cut), path(chromsize)

    output:
    tuple val(meta), path("${meta.id}/summary.*.txt"), optional: true, emit: summary
    tuple val(meta), path("${meta.id}/*.peaks"), optional: true, emit: peak
    path "*.version.txt"          , emit: version

    script:
    def software = "CALL_HIPEAK"
    """
    install_packages.r VGAM MASS
    regression_and_peak_caller.r ${meta.id} $r2peak $r1peak pairs $fasta $mappability $cut $options.args
    echo '1.0.0' > ${software}.version.txt
    """
}
