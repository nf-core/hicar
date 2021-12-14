// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process R1READS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gawk:5.1.0"
    } else {
        container "quay.io/biocontainers/gawk:5.1.0"
    }

    input:
    tuple val(meta), path(pair)

    output:
    tuple val(meta), path("*.bed.gz"), emit: bed
    path "versions.yml"           , emit: versions

    script:
    def software = "awk"
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    gunzip -c $pair | \\
        awk 'BEGIN {OFS="\t"};  /^[^#]/ { print \$2, \$3, \$3+1, "*", "*", \$6}' | \\
        sort -k1,1 -k2,2n | \\
        gzip -nc > ${prefix}.R1.distal.bed.gz

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
    END_VERSIONS
    """
}
