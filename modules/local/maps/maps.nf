process MAPS_MAPS{
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe, stageAs: "long/*"), path(short_bed, stageAs: "short/*"), path(background)
    path make_maps_runfile_source
    val long_bedpe_postfix
    val short_bed_postfix

    output:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe), path(short_bed), path(background), path("${meta.id}_${bin_size}/*"), emit: maps
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
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
        $args

    ## step 2, parse the signals into .xor and .and files, details please refer: doi:10.1371/journal.pcbi.1006982
    ## by default, the sex chromosome will be excluded.
    MAPS.py "${meta.id}_${bin_size}/maps_${meta.id}.maps" $long_bedpe_postfix $short_bed_postfix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
