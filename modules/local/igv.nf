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
    tuple val(name), path(track)
    val species

    output:
    path "{index.html,readme.txt}"   , emit: igv
    path "track.tgz", optional:true  , emit: raw
    path "versions.yml"              , emit: versions

    script:
    def sampleLabel = name.join(' ')
    def tracks      = track.join(' ')
    """
    n=($sampleLabel)
    s=($tracks)
    for i in \$(seq 0 \$(( \${#n[@]} - 1 ))); do
        echo "\${n[\$i]}\t\${s[\$i]}" >> track_files.txt
    done
    create_igv.py track_files.txt $species
    echo "untar the track.tgz and copy all the files into same folder in a web-server." > readme.txt
    tar chvzf track.tgz $track

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS
    """
}
