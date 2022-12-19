process JUICER_ADDNORM {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(hic)
    val juicer_jvm_params
    path juicer_box_jar

    output:
    tuple val(meta), path(hic, includeInputs: true) , emit: hic
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ## add norm just in case there is no normalization data available
    ## the pipeline will using the `-n` parameter when create the .hic if there is only one chromosome

    java ${juicer_jvm_params} \\
        -jar ${juicer_box_jar} \\
        addNorm \\
        --threads $task.cpus \\
        -w ${meta.bin} \\
        $args \\
        $hic

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${juicer_box_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
