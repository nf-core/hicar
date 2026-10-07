process MAPS_FEND {
    tag "$bin_size"
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(bin_size), path(cut)
    path chrom_sizes

    output:
    tuple val(bin_size), path("*.bed")      , emit: bed
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    awk -vOFS="\t" '{print \$3,\$4,\$4,\$3"_"\$1,"0",\$2}' $cut | \\
        bedtools slop $args \\
            -r $bin_size -g $chrom_sizes > \\
            ${cut}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version) | sed -e "s/bedtools v//g")
        MAPS: 1.1.0
    END_VERSIONS
    """
}
