process PAIRSPLOT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bioconductor-chipqc=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-chipqc:1.28.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-chipqc:1.28.0--r41hdfd78af_0"
    }

    input:
    tuple val(meta), path(qc, stageAs: "pairsqc_report/*")

    output:
    tuple val(meta), path("${meta.id}_report/*"), emit: qc
    tuple val(meta), path("${meta.id}_report/*.summary.out"), emit: summary
    tuple val(meta), path("${meta.id}_report/*.distance.vs.proportion.csv"), emit: csv
    path "versions.yml"                        , emit: versions

    script:
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
    mv pairsqc_report ${meta.id}_report
    pairsqcplot.r $enzyme ${meta.id}_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairsqc: "0.2.2"
    END_VERSIONS
    """
}
