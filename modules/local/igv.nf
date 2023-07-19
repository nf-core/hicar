process IGV {
    label 'process_low'
    label 'error_ignore'

    conda "python=3.8"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8' :
        'biocontainers/python:3.8' }"

    input:
    path track
    val species
    val outdir

    output:
    path "{index.html,readme.txt}"   , emit: igv
    path "versions.yml"              , emit: versions

    script:
    """
    create_igv.py $track $species "${outdir.replaceAll('/$', '')}"
    echo "collect in ${track} and copy all the files into relative folder in a web-server." > readme.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS
    """
}
