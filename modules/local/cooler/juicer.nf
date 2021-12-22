process JUICER {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::java-jdk=8.0.112" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1"
    } else {
        container "quay.io/biocontainers/java-jdk:8.0.112--1"
    }

    input:
    tuple val(meta), path(gi)
    path juicer_tools_jar
    path chromsize
    val juicer_jvm_params

    output:
    tuple val(meta), path("*.hic")               , emit: hic
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = task.ext.prefix ? "${meta.id}${task.ext.prefix}" : "${meta.id}"
    def args = task.ext.args ?: ''
    """
    ## thanks https://www.biostars.org/p/360254/
    tail -n +2 $gi | sort -k2,2d -k6,6d > ${gi}.sorted
    # count available chromsomes in the file
    skip_do_norm=\$(awk '{ a[\$2]++; a[\$6]++ } END { if(length(a)==1) print("-n") }' ${gi}.sorted)
    resolutions=(500 1000 2000 5000 10000 20000 50000 100000 250000 500000 1000000 2000000 5000000)
    available_size=(\$(awk '{print \$NF}' $chromsize))
    max_res=5000000
    for i in \${available_size[@]}; do
        if [[ \$i -lt \${max_res} ]]; then
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
    java ${juicer_jvm_params} -jar ${juicer_tools_jar} pre \\
        -r \$res \\
        \${skip_do_norm} \\
        $args --threads $task.cpus ${gi}.sorted ${prefix}.${meta.bin}.hic $chromsize

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(echo \$(java -jar ${juicer_tools_jar} --version 2>&1) | sed 's/^.*Version //; s/Usage.*\$//')
    END_VERSIONS
    """
}
