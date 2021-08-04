// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process GENOME_FILTER {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_1"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_1"
    }

    input:
    path sizes
    path blacklist

    output:
    path "*.bed",         emit: bed
    path "*.version.txt", emit: version

    script:
    def software = "bedtools"
    def file_out = "${sizes.simpleName}.include_regions.bed"
    if (params.blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        echo \$(bedtools --version) | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        echo \$(bedtools --version) | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    }
}
