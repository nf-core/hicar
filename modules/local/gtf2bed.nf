process GTF2BED {
    tag '$gtf'
    label 'process_low'

    conda "bioconda::perl-getopt-long=2.50"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-getopt-long:2.50--pl526_1' :
        'biocontainers/perl-getopt-long:2.50--pl526_1' }"

    input:
    path gtf

    output:
    path "*.bed"                  , emit: bed
    path "versions.yml"           , emit: versions

    script:
    """
    gtf2bed $gtf > ${gtf.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print \$];')
    END_VERSIONS
    """
}
