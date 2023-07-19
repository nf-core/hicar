process PEAKACHU_SCORE {
    label 'process_single'

    conda "bioconda::cooltools=0.5.2"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.2--py39h5371cbf_1' :
        'biocontainers/cooltools:0.5.2--py39h5371cbf_1' }"

    input:
    tuple val(meta), path(cool), path(model)

    output:
    tuple val(meta), val(meta.bin), path("${prefix}.bedpe")    , emit: interactions
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_bin${meta.bin}"
    """
    peakachu_score_genome.py \\
        --resolution ${meta.bin} \\
        --path ${cool} \\
        --model ${model} \\
        --output ${prefix}.bedpe \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peakachu: 2.1
    END_VERSIONS
    """
}
