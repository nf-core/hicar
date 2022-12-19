process HICEXPLORER_HICDETECTLOOPS {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), val(meta.bin), path("${prefix}.bedgraph")   , emit:compartments
    path("versions.yml")                                         , emit:versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    hicDetectLoops \\
        $args \\
        -matrix ${cool} \\
        --outFileName ${prefix}.bedgraph \\
        --threads 1 \\
        --threadsPerChromosome ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicDetectLoops --version 2>&1 | sed 's/hicDetectLoops //')
    END_VERSIONS
    """
}
