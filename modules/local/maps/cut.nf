process MAPS_CUT {
    tag "$bin_size"
    label 'process_high'
    label 'process_long'
    label 'error_retry'

    conda (params.enable_conda ? "conda-forge::biopython=1.70" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_1' :
        'quay.io/biocontainers/biopython:1.70--np112py36_1' }"

    input:
    path fasta
    val bin_size

    output:
    tuple val(bin_size), path('*.cut')        , emit: cut
    path "versions.yml"                       , emit: versions

    script:
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MAPS: 1.1.0
    END_VERSIONS
    """
}
