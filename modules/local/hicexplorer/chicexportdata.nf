process HICEXPLORER_CHICEXPORTDATA {
    tag "${bin_size}"
    label 'process_medium'

    conda "bioconda::hicexplorer=3.7.2 cleanlab=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicexplorer:3.7.2--pyhdfd78af_1' :
        'biocontainers/hicexplorer:3.7.2--pyhdfd78af_1' }"

    input:
    tuple val(bin_size), path(hdf5)

    output:
    tuple val(bin_size), path("${prefix}")      , emit: expdata
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${hdf5.singleName}.tar.gz"
    """
    chicExportData \\
        --file hdf5 \\
        -o $prefix \\
        -t $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hicexplorer: \$(chicExportData --version 2>&1 | sed 's/chicExportData //')
    END_VERSIONS
    """
}
