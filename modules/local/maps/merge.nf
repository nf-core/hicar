// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_MERGE {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    tuple val(bin_size), path(cut), path(mappability)
    path merge_map_py_source

    output:
    tuple val(bin_size), path("${cut.getSimpleName()}")    , emit: map
    path "versions.yml"                                    , emit: versions

    script:
    """
    python $merge_map_py_source \\
        -c $cut \\
        -m $mappability \\
        -o tmp.map
    awk "\\\$7>${options.args}" tmp.map > ${cut.getSimpleName()}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        MAPS: 1.1.0
    END_VERSIONS
    """
}
