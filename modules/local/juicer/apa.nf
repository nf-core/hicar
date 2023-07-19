process JUICER_APA {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(hic), path(loops)
    path juicer_tools_jar

    output:
    tuple val(meta), path("$prefix/*")           , emit: results
    tuple val(meta), path("$prefix/**.png")      , emit: plot
    path "versions.yml"                          , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def avail_mem = 4
    if (!task.memory) {
        log.info 'Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    java -Xms512m -Xmx${avail_mem}g -jar ${juicer_tools_jar} apa \\
        $args \\
        --threads $task.cpus \\
        $hic $loops $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${juicer_tools_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
