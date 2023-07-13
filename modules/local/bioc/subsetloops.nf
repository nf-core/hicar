process BIOC_SUBSETLOOPS {
    tag "$mergedloops"
    label 'process_single'

    conda "bioconda::bioconductor-rtracklayer=1.50.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' :
        'biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2' }"

    input:
    path peaks
    tuple val(bin_size), path(mergedloops)

    output:
    tuple val(bin_size), path("sub_loop/${prefix}*")    , emit: loops
    tuple val(bin_size), path("sub_peak/${prefix}*")    , emit: bed
    path "versions.yml"                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'subset_'
    """
    subsetloops.r \\
        --task_process ${task.process} \\
        --prefix $prefix \\
        --peaks $peaks \\
        $mergedloops

    # versions.yml files will be created in the rscripts
    """
}
