process RE_CUTSITE {
    label 'process_low'

    conda "conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.78' :
        'quay.io/biocontainers/biopython:1.78' }"

    input:
    val enzyme

    output:
    stdout emit: site
    path "versions.yml"              , emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version) | sed 's/Python //')
    END_VERSIONS

    restriction_enzyme_cutsite.py $enzyme
    """
}
