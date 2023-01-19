process JUICER_EIGENVECTOR {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(hic)
    val resolution
    tuple path(juicer_box_jar), path(chromsizes)

    output:
    tuple val(meta), path("$prefix/*")           , emit: compartments
    path "versions.yml"                          , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}_${resolution}"
    def args = task.ext.args ?: ''
    norm_method = (args.contains('NONE')) ? 'NONE' :
        (args.contains('VC_SQRT')) ? 'VC_SQRT' :
        (args.contains('VC')) ? 'VC' :
        (args.contains('KR')) ? 'KR' : 'NONE'
    data_slot = (args.contains('BP')) ? 'BP' :
        (args.contains('FRAG')) ? 'FRAG' : 'BP'
    args = args.tokenize()
    args.removeIf { it.contains('NONE') }
    args.removeIf { it.contains('VC_SQRT') }
    args.removeIf { it.contains('VC') }
    args.removeIf { it.contains('KR') }
    args.removeIf { it.contains('BP') }
    args.removeIf { it.contains('FRAG') }
    args = args.join(' ')
    def avail_mem = 4
    if (!task.memory) {
        log.info 'Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    ## add norm just in case there is no normalization data available
    java -Xms512m -Xmx${avail_mem}g \\
        -jar ${juicer_box_jar} \\
        addNorm \\
        --threads $task.cpus \\
        -w ${meta.bin} \\
        -k $norm_method \\
        $hic
    mkdir -p ${prefix}
    while read -r chrom; do
        chrom=\${chrom%\$'\\t'*}
        java -Xms512m -Xmx${avail_mem}g \\
            -jar ${juicer_box_jar} \\
            eigenvector \\
            --threads $task.cpus \\
            $args \\
            $norm_method \\
            $hic \\
            \$chrom \\
            $data_slot \\
            $resolution \\
            ${prefix}/${prefix}_\${chrom}.eigen.txt
    done < $chromsizes


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${juicer_box_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
