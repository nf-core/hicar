// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IGV {
    label 'process_low'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "python=3.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8"
    } else {
        container "quay.io/biocontainers/python:3.8"
    }

    input:
    path track
    val species

    output:
    path "{index.html,readme.txt}"   , emit: igv
    path "versions.yml"              , emit: versions

    script:
    """
    create_igv.py $track $species
    echo "collect in ${track} and copy all the files into relative folder in a web-server. Then put the index.html file into the root folder." > readme.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS
    """
}
