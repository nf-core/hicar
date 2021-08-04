// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GTF2BED {
    tag '$gtf'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) },
        enabled: options.publish

    conda (params.conda ? "conda-forge::perl=5.26.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sed:4.2.3.dev0--0"
    } else {
        container "quay.io/biocontainers/sed:4.2.3.dev0--0"
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
