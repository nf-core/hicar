// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_MAPS{
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe, stageAs: "long/*"), path(short_bed, stageAs: "short/*"), path(background)
    path make_maps_runfile_source

    output:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe), path(short_bed), path(background), path("${meta.id}_${bin_size}/*"), emit: maps
    path "versions.yml"          , emit: versions

    script:
    def software  = "MAPS"
    """
    ## 2 steps
    ## step 1, prepare the config file for MAPS. The file will be used for multiple steps
    mkdir -p "${meta.id}_${bin_size}"
    python ${make_maps_runfile_source} \\
        "${meta.id}" \\
        "${meta.id}_${bin_size}/" \\
        $macs2 \\
        $background \\
        "long/" \\
        "short/" \\
        $bin_size \\
        0 \\
        "${meta.id}_${bin_size}/" \\
        ${options.args}
    ## step 2, parse the signals into .xor and .and files, details please refer: doi:10.1371/journal.pcbi.1006982
    ## by default, the sex chromosome will be excluded.
    MAPS.py "${meta.id}_${bin_size}/maps_${meta.id}.maps"

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        MAPS: 1.1.0
    END_VERSIONS
    """
}
