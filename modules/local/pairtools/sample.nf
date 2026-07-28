process PAIRTOOLS_SAMPLE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pairtools=0.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:0.3.0--py37hb9c2fc3_5' :
        'biocontainers/pairtools:0.3.0--py37hb9c2fc3_5' }"

    input:
    tuple val(meta), path(input), val(fraction)

    output:
    tuple val(meta), path("*.pairsam.gz")  , emit: pairsam
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_subsample"
    """
    pairtools \\
        sample \\
        $args \\
        --nproc-in $task.cpus \\
        --nproc-out $task.cpus \\
        -o ${prefix}.pairsam.gz \\
        $fraction \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pairtools: \$(pairtools --version 2>&1 | sed 's/pairtools.*version //')
    END_VERSIONS
    """
}
