// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_PEAK {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
    path peak

    output:
    path "merged_peak.bed"    , emit: peak
    path "versions.yml"       , emit: versions

    script:
    def software = "bedtools"
    """
    cat *.narrowPeak | \\
        cut -f1-3 | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge $options.args \\
            -i stdin > merged_peak.bed

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bedtools --version) | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
