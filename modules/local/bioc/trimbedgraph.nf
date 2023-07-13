process BEDGRAPH_TRIM {
    tag "$bedgraph"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    tuple val(meta), path(bedgraph)
    path  sizes

    output:
    tuple val(meta), path("*.trimmed.bedgraph")         , emit: bedgraph
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: 'bedGraph'
    """
    trimbedgraph.r \\
        --format $args \\
        --chrom_size $sizes \\
        $bedgraph

    # *.version.txt files will be created in the rscripts
    """
}
