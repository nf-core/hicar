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
    path juicer_box_jar

    output:
    tuple val(meta), path(hic, includeInputs: true) , emit: hic
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 4
    if (!task.memory) {
        log.info 'Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    ## add norm just in case there is no normalization data available
    ## the pipeline will using the `-n` parameter when create the .hic if there is only one chromosome

    java -Xms512m -Xmx${avail_mem}g \\
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
