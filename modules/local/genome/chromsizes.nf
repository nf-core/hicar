process CHROMSIZES {
    tag "$fasta"
    label 'process_low'

    conda "bioconda::samtools=1.12"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0' :
        'biocontainers/samtools:1.12--hd5e65b6_0' }"

    input:
    path fasta

    output:
    path '*.sizes'      , emit: sizes
    path '*.fai'        , emit: fai
    path "versions.yml" , emit: versions

    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai | sort -k 1,1 > ${fasta}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
