process COOLTOOLS_COMPARTMENTS {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cooltools=0.5.1 bioconda::ucsc-bedgraphtobigwig=377" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' :
        'quay.io/biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' }"

    input:
    tuple val(meta), path(mcool)
    val resolution
    path fasta
    path chromsizes

    output:
    tuple val(meta), path("*compartments*")         , emit: results
    tuple val(meta), path('*.bw')                   , emit: compartments
    path("versions.yml")                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cooltools genome binnify --all-names ${chromsizes} ${resolution} > genome_bins.txt
    cooltools genome gc genome_bins.txt ${fasta} > genome_gc.txt
    cooltools eigs-cis \\
        $args \\
        --phasing-track genome_gc.txt \\
        -o ${prefix}_compartments \\
        ${mcool}::resolutions/${resolution}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | sed 's/cooletools, version //')
    END_VERSIONS
    """
}
