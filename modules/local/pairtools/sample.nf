process PAIRTOOLS_SAMPLE {
    tag "$meta.id"
    label 'process_high'

    // Pinning numpy to 1.23 until https://github.com/open2c/pairtools/issues/170 is resolved
    // Not an issue with the biocontainers because they were built prior to numpy 1.24
    conda "bioconda::pairtools=1.0.2 conda-forge::numpy=1.23"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pairtools:1.0.2--py39h2a9f597_0' :
        'quay.io/biocontainers/pairtools:1.0.2--py39h2a9f597_0' }"

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
        pairtools: \$(pairtools --version 2>&1 | tr '\\n' ',' | sed 's/.*pairtools.*version //' | sed 's/,\$/\\n/')
    END_VERSIONS
    """
}
