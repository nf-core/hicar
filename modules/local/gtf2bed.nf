// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GTF2BED {
    tag '$gtf'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

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
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    gtf2bed $gtf > ${gtf.baseName}.bed

    echo \$(perl -e 'print \$];') > perl.version.txt
    """
}
