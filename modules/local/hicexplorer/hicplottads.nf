process HICEXPLORER_HICPLOTTADS {
    tag "${meta.id}"
    label 'process_high'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::hicexplorer=3.7.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'quay.io/biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(cool), path(tads)
    val chrom_sizes

    output:
    tuple val(meta), path("*.png")              , emit:pngs
    tuple val(meta), path("*.ini")              , emit:ini
    path("versions.yml")                        , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.bin}"
    """
    #step 1, create the hicar_track.ini
    create_hicar_track_ini.py \\
        --matrix $cool \\
        --tads $tads \\
        --out ${prefix}_hicar_track.ini
    while read chr seqlength; do
        #step 2, plot the TADs
        hicPlotTADs \\
            ${args} \\
            --tracks ${prefix}_hicar_track.ini \\
            -o ${prefix}_track_\${chr}.png \\
            --region \${chr}:1-\$seqlength

    done < ${chrom_sizes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(hicPlotTADs --version 2>&1 | sed 's/hicPlotTADs //')
    END_VERSIONS
    """
}
