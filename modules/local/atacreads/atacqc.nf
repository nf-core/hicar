// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process ATACQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bioconductor-atacseqqc=1.16.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-atacseqqc:1.16.0--r41hdfd78af_0"
    }

    input:
    path peaks
    path beds
    path gtf

    output:
    path "*.csv"                   , emit: stats
    path  "*.version.txt"          , emit: version

    script:
    def software = "ATACseqQC"
    """
    atacqc.r $gtf

    # *.version.txt files will be created in the rscripts
    """
}
