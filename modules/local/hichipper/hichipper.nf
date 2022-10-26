process HICHIPPER_HICHIPPER {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::hichipper=0.7.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hichipper:0.7.7--py_1' :
        'quay.io/biocontainers/hichipper:0.7.7--py_1' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), path("*${prefix}.h5")                    , emit:interactions
    path("versions.yml")                                      , emit:versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    ## step1 create yaml file1
    ## step2 run programs
    hichipper \\
        -m $cool \\
        $args \\
        --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hichipper: \$(hichipper --version 2>&1 | sed 's/hichipper //')
    END_VERSIONS
    """
}
