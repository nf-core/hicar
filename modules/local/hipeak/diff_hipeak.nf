// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process DIFF_HIPEAK {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-diffhic=1.24.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-diffhic:1.24.0--r41h399db7b_0 "
    } else {
        container "quay.io/biocontainers/bioconductor-diffhic:1.24.0--r41h399db7b_0"
    }

    input:
    path peaks, stageAs: "peaks/*"
    path summary, stageAs: "summary/*"
    path distalpair, stageAs: "pairs/*"

    output:
    path "${prefix}/*"                        , emit: diff
    path "${prefix}/*.qc.json", optional: true, emit: stats
    path "*.version.txt"                      , emit: version

    script:
    prefix   = options.suffix ? "${options.suffix}" : "diffhicar"
    """
    diff_hipeak.r $prefix \\
        $options.args
    # *.version.txt files will be created in the rscripts
    """
}
