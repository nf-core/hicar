// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process ENSEMBL_UCSC_CONVERT {
    tag "$fname"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bioconductor-rtracklayer=1.50.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    } else {
        container "quay.io/biocontainers/bioconductor-rtracklayer:1.50.0--r40h7f5ccec_2"
    }

    input:
    tuple val(bin_size), path(fname)

    output:
    tuple val(bin_size), path("{UCSC,ensembl}.${fname}"), emit: tab
    path "versions.yml"                                 , emit: versions

    script:
    """
    seqlevels_convert.r \\
        $options.args \\
        $fname

    # *.version.txt files will be created in the rscripts
    echo "${getProcessName(task.process)}:" > versions.yml
    for i in \$(ls *.version.txt); do
    echo "    \${i%.version.txt}: \$(<\$i)" >> versions.yml
    done
    """
}
