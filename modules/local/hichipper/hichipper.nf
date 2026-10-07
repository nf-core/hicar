process HICHIPPER_HICHIPPER {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::hichipper=0.7.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hichipper:0.7.7--py_1' :
        'biocontainers/hichipper:0.7.7--py_1' }"

    input:
    tuple val(meta), path(hicpro), path(peak), path(yaml)
    path resfrags

    output:
    tuple val(meta), path("*${prefix}")                       , emit:interactions
    path("versions.yml")                                      , emit:versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    ## step2 run programs
    hichipper \\
        --out ${prefix} \\
        $yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hichipper: \$(hichipper --version 2>&1 | sed 's/hichipper //')
    END_VERSIONS
    """
}
