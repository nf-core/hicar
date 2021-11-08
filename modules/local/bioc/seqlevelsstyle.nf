// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process SEQLEVELS_STYLE {
    tag "$bed"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::bioconductor-genomeinfodb=1.26.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0"
    } else {
        container "quay.io/biocontainers/bioconductor-genomeinfodb:1.26.4--r40hdfd78af_0"
    }

    input:
    path bed

    output:
    stdout emit: seqlevels_style
    path "versions.yml"           , emit: versions

    script:
    def software = "GenomeInfoDb"

    """
    seqlevelsstyle.r $bed > r.log.txt 2>&1

    # *.version.txt files will be created in the rscripts
    echo "${getProcessName(task.process)}:" > versions.yml
    for i in \$(ls *.version.txt); do
    echo "    \${i%.version.txt}: \$(<\$i)" >> versions.yml
    done
    """
}
