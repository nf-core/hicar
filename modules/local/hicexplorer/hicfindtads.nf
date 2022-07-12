process HICEXPLORER_HICFINDTADS {
    tag "${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::hicexplorer=3.7.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("*hicfindtads*")      , emit:results
    tuple val(meta), path("*domains.bed")       , emit:domains
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    def bin_size = meta.bin.toInteger()
    def minDepth = bin_size * 3
    def maxDepth = bin_size * 10
    """
    hicFindTADs \\
        ${args} \\
        --matrix ${cool} \\
        --outPrefix ${prefix}_hicfindtads \\
        --minDepth $minDepth \\
        --maxDepth $maxDepth \\
        --step $resolution \\
        --numberOfProcessors ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicFindTADs --version 2>&1 | sed 's/hicFindTADs //')
    END_VERSIONS
    """
}
