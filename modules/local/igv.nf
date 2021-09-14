// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

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
    tuple val(name), path(track)
    val species

    output:
    path "{index.html,readme.txt}"   , emit: igv
    path "track.tgz", optional:true  , emit: raw
    path "*.version.txt"             , emit: version

    script:
    def sampleLabel = name.join(' ')
    def tracks      = track.join(' ')
    """
    n=($sampleLabel)
    s=($tracks)
    paste <(printf '%s\n' "\${n[@]}") <(printf '%s\n' "\${s[@]}") > track_files.txt
    create_igv.py track_files.txt $species
    echo "untar the track.tgz and copy all the files into same folder in a web-server." > readme.txt
    tar chvzf track.tgz $track

    echo \$(python --version) | sed 's/Python //'> python.version.txt
    """
}
