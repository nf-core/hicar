// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_MERGE {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    tuple val(bin_size), path(cut), path(mappability)

    output:
    tuple val(bin_size), path("${cut.getSimpleName()}")    , emit: map
    path "*.version.txt"                                   , emit: version

    script:
    def software = "MAPS"
    """
    merge_map.py \\
        -c $cut \\
        -m $mappability \\
        -o tmp.map
    awk '\$7>0.5' tmp.map > ${cut.getSimpleName()}

    echo '1.1.0' > ${software}.version.txt
    """
}
