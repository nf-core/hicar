process HICEXPLORER_HICAGGREGATECONTACTS {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::hicexplorer=3.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool)
    path anchor

    output:
    tuple val(meta), path("*apa*")              , emit:results
    tuple val(meta), path("*.png")              , emit:png
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    hicAggregateContacts \\
        ${args} \\
        --matrix ${cool} \\
        --BED ${anchor} \\
        --outFileName ${prefix}_apa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicFindTADs --version 2>&1 | sed 's/hicFindTADs //')
    END_VERSIONS
    """
}
