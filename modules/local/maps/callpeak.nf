process MAPS_CALLPEAK {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "bioconda::bioconductor-monocle=2.20.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-monocle:2.20.0--r41h399db7b_0' :
        'biocontainers/bioconductor-monocle:2.20.0--r41h399db7b_0' }"

    input:
    tuple val(meta), val(bin_size), path(macs2), path(long_bedpe, stageAs: "long/*"), path(short_bed, stageAs: "short/*"), path(background), path(maps, stageAs: "maps_out/*")

    output:
    tuple val(meta), val(bin_size), path("${meta.id}_${bin_size}/summary.*.txt"), optional: true, emit: summary
    tuple val(meta), val(bin_size), path("${meta.id}_${bin_size}/*.peaks"), optional: true, emit: peak
    tuple val(meta), val(bin_size), path("${meta.id}_${bin_size}/reg_raw.*.MAPS2_*"), optional: true, emit: signal
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mv maps_out ${meta.id}_${bin_size}
    ## arguments:
    ## INFDIR - dir with reg files, output folder of maps/maps.nf::MAP_MAPS
    ## SET - dataset name, output of maps/maps.nf::MAP_MAPS. It is the postfix of the file name of reg files, eg: reg_raw.chr22.WT.5k.xor, SET should be WT.5k.
    ## RESOLUTION - resolution (for example 5000 or 10000). In this pipeline, it is the bin_size.
    ## COUNT_CUTOFF - count cutoff, default 12
    ## RATIO_CUTOFF - ratio cutoff, default 2.0
    ## FDR - -log10(fdr) cutoff, default 2
    ## FILTER - file containing bins that need to be filtered out. Format: two columns "chrom", "bin". "chrom" contains 'chr1','chr2',.. "bin" is bin label
    ## regresison_type - pospoisson for positive poisson regression, negbinom for negative binomial. default is pospoisson
    MAPS_regression_and_peak_caller.r "${meta.id}_${bin_size}/" ${meta.id} $bin_size $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    for i in \$(ls *.version.txt); do
    echo "    \${i%.version.txt}: \$(<\$i)" >> versions.yml
    done
    """
}
