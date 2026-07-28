process JUICER_PRE {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(gi)
    path hic_tools_jar
    path chromsize

    output:
    tuple val(meta), path("*.hic")               , emit: hic
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def avail_mem = 4
    if (!task.memory) {
        log.info 'Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    ## thanks https://www.biostars.org/p/360254/
    resolutions=(500 1000 2000 5000 10000 20000 50000 100000 250000 500000 1000000 2000000 5000000)
    ## clean the chromosomes, otherwise, the patch chromosomes size may be very small
    ## skip this line if it does not fit your model system
    grep -v [_M] $chromsize > ${chromsize}tmp
    grep -v EBV ${chromsize}tmp > ${chromsize}_fil.txt
    rm ${chromsize}tmp
    nl=\$(wc -l < ${chromsize}_fil.txt)
    if [ \$nl -eq 0 ]; then
        cp $chromsize ${chromsize}_fil.txt
    fi
    available_size=(\$(awk '{print \$NF}' ${chromsize}_fil.txt))
    max_res=5000000
    for i in \${available_size[@]}; do
        if [ \$i -lt \${max_res} ]; then
            max_res=\$i
        fi
    done
    res=()
    function join_by { local IFS="\$1"; shift; echo "\$*"; }
    for i in \${resolutions[@]}
    do
        if [ \$i -ge $meta.bin ] && [ \$i -lt \${max_res} ]; then
            res+=(\$i)
        fi
    done

    res=\$(join_by , \${res[@]})

    ## sort and filter the input pairs
    seqlev=(\$(awk '{print \$1}' ${chromsize}_fil.txt))
    seqlev=\$(join_by '|' \${seqlev[@]})
    tail -n +2 $gi | \\
        awk -v seqlev=\$seqlev 'match(\$2, seqlev) && match(\$6, seqlev) {print}' | \\
        sort -k2,2d -k6,6d -k3,3n -k7,7n > ${gi}.sorted
    # count available chromsomes in the file
    skip_do_norm=\$(awk '{ a[\$2]++; a[\$6]++ } END { if(length(a)==1) print("-n") }' ${gi}.sorted)

    java -Xms512m -Xmx${avail_mem}g \\
        -jar ${hic_tools_jar} pre \\
        -r \$res \\
        \${skip_do_norm} \\
        $args \\
        --threads $task.cpus \\
        ${gi}.sorted ${prefix}.${meta.bin}.hic ${chromsize}_fil.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${hic_tools_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
