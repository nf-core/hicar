// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process PAIRSPLOT {
    tag "$meta.id"
    label 'process_low'
    errorStrategy { (task.exitStatus in 137..140 && task.attempt <= 3)  ? 'retry' : 'ignore' }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::r-nozzle.r1=1.1_1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-nozzle.r1:1.1_1--r3.3.2_0"
    } else {
        container "quay.io/biocontainers/r-nozzle.r1:1.1_1--r3.3.2_0"
    }

    input:
    tuple val(meta), path(qc, stageAs: "pairsqc_report/*")

    output:
    tuple val(meta), path("${meta.id}_report/*"), emit: qc
    path "*.version.txt"                        , emit: version

    script:
    def software = "pairsqc"
    def RE_cutsite = [
        "mboi": 4,
        "ncoi": 6,
        "dpnii": 4,
        "bglii": 6,
        "hindiii": 6,
        "cviqi": 4,
        "arima": 4,
        "mnase": 4]
    def enzyme = RE_cutsite[params.enzyme.toLowerCase()]?:4
    """
    install_packages.r SooLee/plotosaurus stringr Nozzle.R1
    mv pairsqc_report ${meta.id}_report
    pairsqcplot.r $enzyme ${meta.id}_report

    echo "0.2.2" > ${software}.version.txt
    """
}
