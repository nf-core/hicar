// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_FEATURE {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "bioconda::samtools=1.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.10--h9402c20_2"
    } else {
        container "quay.io/biocontainers/samtools:1.10--h9402c20_2"
    }

    input:
    tuple val(bin_size), path(map)
    path chrom_sizes

    output:
    tuple val(bin_size), path("*_el.txt")  , emit: bin_feature
    path "*.version.txt"                   , emit: version

    script:
    def software = "MAPS"
    """
    feature_frag2bin.py \\
        -i $map \\
        -o F_GC_M_${map.getSimpleName()}_${bin_size}_el.txt \\
        -b $bin_size \\
        -g $chrom_sizes

    echo '1.1.0' > ${software}.version.txt
    """
}
