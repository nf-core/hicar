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
    java -Xms512m -Xmx2048m -jar $projectDir/bin/juicer_tools_1.22.01.jar pre \
        $options.args --threads $task.cpus ${gi}.sorted ${prefix}.hic $chromsize

    echo 1.22.01 > ${software}.version.txt
    """
}
