// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_REFORMAT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::r-data.table=1.12.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-data.table:1.12.2"
    } else {
        container "quay.io/biocontainers/r-data.table:1.12.2"
    }

    input:
    tuple val(meta), val(bin_size), path(peak)

    output:
    tuple val(meta), val(bin_size), path("*.sig3Dinteractions.bedpe"), emit: bedpe
    path "versions.yml"           , emit: versions

    script:
    """
    MAPS_peak_formatting.r $bin_size ${peak}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        MAPS: 1.1.0
    END_VERSIONS
    for i in \$(ls *.version.txt); do
    echo "    \${i%.version.txt}: \$(<\$i)" >> versions.yml
    done
    """
}
