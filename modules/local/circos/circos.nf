process CIRCOS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::circos=0.69.8" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/circos:0.69.8--hdfd78af_1' :
        'quay.io/biocontainers/circos:0.69.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(data)
    path configfile

    output:
    path "*.png"                  , emit: circos
    path "versions.yml"           , emit: versions

    script:
    """
    circos
    mv circos.png ${meta.id}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circos: \$(echo \$(circos -v 2>&1) | sed 's/circos.*v //; s/ .*\$//')
    END_VERSIONS
    """
}
