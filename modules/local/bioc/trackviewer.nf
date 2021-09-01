// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BIOC_TRACKVIEWER {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-trackviewer=1.28.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-trackviewer:1.28.0--r41h399db7b_0"
    }

    when:
    !workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container

    input:
    tuple val(bin_size), path(events), path(mcools)
    path raw_pairs // .unselected.pairs.gz of samfrag
    path gtf
    path chrom_sizes
    path restrict

    output:
    path "${prefix}/*"            , emit: v4c
    path "*.version.txt"          , emit: version

    script:
    prefix   = options.suffix ? "${options.suffix}${bin_size}" : "diffhic_bin${bin_size}"
    """
    install_packages.r trackViewer optparse
    trackviewer.r -g "${gtf}" -s "${chrom_sizes}" -x "${restrict}" -r $bin_size -o "${prefix}" $options.args
    """
}
