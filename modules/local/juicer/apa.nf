process JUICER_APA {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::java-jdk=8.0.112" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(hic)
    path loops
    path juicer_box_jar
    val juicer_jvm_params

    output:
    tuple val(meta), path("$prefix/*")           , emit: results
    tuple val(meta), path("$prefix/*.png")       , emit: png
    path "versions.yml"                          , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    java ${juicer_jvm_params} -jar ${juicer_box_jar} apa \\
        $args \\
        --threads $task.cpus \\
        $hic $loops $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${juicer_tools_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
