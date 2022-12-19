process MAPS_FEATURE {
    tag "$bin_size"
    label 'process_low'

    conda "pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(bin_size), path(map)
    path chrom_sizes
    path feature_frag2bin_source

    output:
    tuple val(bin_size), path("*_el.txt")  , emit: bin_feature
    path "versions.yml"                    , emit: versions

    script:
    """
    python ${feature_frag2bin_source} \\
        -i $map \\
        -o F_GC_M_${map.getSimpleName()}_${bin_size}_el.txt \\
        -b $bin_size \\
        -g $chrom_sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
