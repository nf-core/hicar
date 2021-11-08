// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_FEND {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1"
    }

    input:
    tuple val(bin_size), path(cut)
    path chrom_sizes

    output:
    tuple val(bin_size), path("*.bed")      , emit: bed
    path "versions.yml"                     , emit: versions

    script:
    """
    awk -vOFS="\t" '{print \$3,\$4,\$4,\$3"_"\$1,"0",\$2}' $cut | \\
    bedtools slop $options.args \\
        -r $bin_size -g $chrom_sizes > \\
        ${cut}.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        bedtools: \$(echo \$(bedtools --version) | sed -e "s/bedtools v//g")
        MAPS: 1.1.0
    END_VERSIONS
    """
}
