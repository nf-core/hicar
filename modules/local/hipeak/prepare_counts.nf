// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PREPARE_COUNTS {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(r2peak, stageAs:"R2peak/*"), path(r1peak, stageAs: "R1peak/*"), path(distalpair, stageAs: "pairs/*"), path(mappability), path(fasta), path(cut)

    output:
    tuple val(meta), path("*.csv"), emit: counts
    path "*.version.txt"          , emit: version

    script:
    def software = "CALL_HIPEAK"
    """
    install_packages.r rtracklayer InteractionSet Biostrings Rsamtools

    prepare_counts.r --r2peak $r2peak --r1peak $r1peak \\
        --pairs pairs \\
        --fasta $fasta \\
        --mappability $mappability \\
        --restrict $cut \\
        --output counts.${meta.id}.csv \\
        $options.args
    """
}
