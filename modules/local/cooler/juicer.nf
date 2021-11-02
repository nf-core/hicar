// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process JUICER {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::java-jdk=8.0.112" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1"
    } else {
        container "quay.io/biocontainers/java-jdk:8.0.112--1"
    }

    input:
    tuple val(meta), path(gi)
    path chromsize

    output:
    tuple val(meta), path("*.hic")               , emit: hic
    path  "*.version.txt"                        , emit: version

    script:
    def software = "Juicer_tools"
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ## thanks https://www.biostars.org/p/360254/
    sort -k2,2d -k6,6d $gi > ${gi}.sorted
    resolutions=(500 1000 2000 5000 10000 20000 50000 100000 250000 500000 1000000 2000000 5000000)
    res=()
    function join_by { local IFS="\$1"; shift; echo "\$*"; }
    for i in \${resolutions[@]}
    do
        if [[ \$i -ge $meta.bin ]]; then
            res+=(\$i)
        fi
    done
    res=\$(join_by , \${res[@]})
    java ${params.juicer_jvm} -jar $projectDir/bin/juicer_tools_1.22.01.jar pre \
        -r \$res \
        $options.args --threads $task.cpus ${gi}.sorted ${prefix}.${meta.bin}.hic $chromsize

    echo 1.22.01 > ${software}.version.txt
    """
}
