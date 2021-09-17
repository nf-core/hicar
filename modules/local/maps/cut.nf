// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MAPS_CUT {
    tag "$bin_size"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::biopython=1.70" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_1"
    } else {
        container "quay.io/biocontainers/biopython:1.70--np112py36_1"
    }

    input:
    path fasta
    val bin_size

    output:
    tuple val(bin_size), path('*.cut')        , emit: cut
    path "*.version.txt", emit: version

    script:
    def software = 'MAPS'
    def RE_cutsite = [
        "mboi": [site:"GATC", pos:"0"],
        "ncoi": [site:"CCATGG", pos:"1"],
        "dpnii": [site:"GATC", pos:"0"],
        "bglii": [site:"AGATCT", pos:"1"],
        "hindiii": [site:"AAGCTT", pos:'1'],
        "cviqi": [site:"GTAC", pos:"1"],
        "arima": [site:"GATC,GA.TC", pos:"0,1"],
        "mnase": [site:"mnase", pos:"none"]]
    def enzyme = RE_cutsite[params.enzyme.toLowerCase()]?:[site:params.enzyme.replaceAll("^", ""), pos:params.enzyme.indexOf('^')]
    """
    restriction_cut_multipleenzyme.py \\
        -f ${fasta} \\
        -s ${enzyme.site} \\
        -p ${enzyme.pos} \\
        -b ${bin_size} \\
        -o ${bin_size}_${params.enzyme}.cut \\
        -c $task.cpus

    echo '1.1.0' > ${software}.version.txt
    """
}
