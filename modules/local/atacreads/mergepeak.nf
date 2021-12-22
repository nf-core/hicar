process MERGE_PEAK {
    label 'process_medium'

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
    def args = task.ext.args ?: ''
    """
    cat *.narrowPeak | \\
        cut -f1-3 | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge $args \\
            -i stdin > merged_peak.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version) | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
