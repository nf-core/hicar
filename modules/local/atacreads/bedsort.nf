// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process BEDFILES_SORT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::coreutils=8.31" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0 "
    } else {
        container "quay.io/biocontainers/coreutils:8.31--h14c3975_0 "
    }

    input:
    tuple val(meta), path(intervals)
    val   extension

    output:
    tuple val(meta), path("*.${extension}"), emit: sorted
    path  "versions.yml"                   , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def buffer   = task.memory.toGiga().intdiv(2)
    """
    ## ref: https://www.biostars.org/p/66927/
    LC_ALL=C sort \\
        --parallel=$task.cpus \\
        --buffer-size=${buffer}G \\
        -k1,1 -k2,2n \\
        $intervals \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(sort --version | tr '\\n' ' ' | sed -e "s/^[^0-9]*//; s/ Copyright.*\$//")
    END_VERSIONS
    """
}
