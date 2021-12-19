process BIOC_ENRICH {
    tag "$bin_size"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::bioconductor-clusterprofiler=3.18.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-clusterprofiler:3.18.1--r40hdfd78af_0"
    }

    input:
    tuple val(bin_size), path(diff)
    val ucscname

    output:
    tuple val(bin_size), path("${prefix}/enrichment/*"), emit: enrichment
    path "versions.yml"                                , emit: versions

    script:
    prefix   = task.ext.prefix ? "${task.ext.prefix}${bin_size}" : "diffhic_bin${bin_size}"
    """
    enrich.r -s ${ucscname} -o "${prefix}/enrichment" $options.args

    # *.version.txt files will be created in the rscripts
    echo "${task.process}:" > versions.yml
    for i in \$(ls *.version.txt); do
    echo "    \${i%.version.txt}: \$(<\$i)" >> versions.yml
    done
    """
}
