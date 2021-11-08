// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process CIRCOS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::circos=0.69.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/circos:0.69.8--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/circos:0.69.8--hdfd78af_1"
    }

    input:
    tuple val(meta), path(data), path(configfile)

    output:
    path "*.png"                  , emit: circos
    path "versions.yml"           , emit: versions

    script:
    def software = getSoftwareName(task.process)
    """
    circos
    mv circos.png ${meta.id}.png

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(circos -v 2>&1) | sed 's/circos.*v //; s/ .*\$//')
    END_VERSIONS
    """
}
