process GTF2BED {
    tag '$gtf'
    label 'process_low'

    conda (params.enable_conda ? "bioconda::perl-getopt-long=2.50" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl-getopt-long:2.50--pl526_1"
    } else {
        container "quay.io/biocontainers/perl-getopt-long:2.50--pl526_1"
    }

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
