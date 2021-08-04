// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_REFORMAT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::r-vgam=1.0_2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-vgam:1.0_2--r3.3.2_0"
    } else {
        container "quay.io/biocontainers/r-vgam:1.0_2--r3.3.2_0"
    }

    input:
    tuple val(meta), val(bin_size), path(peak)

    output:
    tuple val(meta), val(bin_size), path("*.sig3Dinteractions.bedpe"), emit: bedpe
    path "*.version.txt"          , emit: version

    script:
    def software = "MAPS"
    """
    install_packages.r data.table
    MAPS_peak_formatting.r ${peak} $bin_size

    echo '1.1.0' > ${software}.version.txt
    """
}
