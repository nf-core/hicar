process PEAKACHU_SCORE {
    label 'process_single'

    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(cool), path(model)

    output:
    tuple val(meta), val(meta.bin), path("${prefix}/*")    , emit: interactions
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_bin${meta.bin}"
    """
    python peakachu_score_genome.py \\
        --resolution ${meta.bin} \\
        --path ${cool} \\
        --model ${model} \\
        --output $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peakachu: 2.1
    END_VERSIONS
    """
}
