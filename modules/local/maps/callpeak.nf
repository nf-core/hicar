// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_CALLPEAK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::r-vgam=1.0_2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-vgam:1.0_2--r3.3.2_0"
    } else {
        container "quay.io/biocontainers/r-vgam:1.0_2--r3.3.2_0"
    }

    input:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe, stageAs: "long/*"), path(short_bed, stageAs: "short/*"), path(background), path(maps, stageAs: "maps_out/*")

    output:
    tuple val(meta), val(bin_size), path("${meta.id}_${bin_size}/summary.*.txt"), optional: true, emit: summary
    tuple val(meta), val(bin_size), path("${meta.id}_${bin_size}/*.peaks"), optional: true, emit: peak
    path "*.version.txt"          , emit: version

    script:
    def software = "MAPS"
    def cutoff   = options.args  ?:"12 2.0"
    def sex_chr  = options.args2 ?:""
    def model    = options.args3 ?:"pospoisson"
    def filter   = "None"
    """
    resolution=\$(bc<<<"$bin_size/1000")
    install_packages.r VGAM
    mv maps_out ${meta.id}_${bin_size}
    MAPS_regression_and_peak_caller.r "${meta.id}_${bin_size}/" ${meta.id}.\${resolution}k $bin_size $params.autosomal$sex_chr $filter "${meta.id}_${bin_size}/" $cutoff $model
    echo '1.1.0' > ${software}.version.txt
    """
}
