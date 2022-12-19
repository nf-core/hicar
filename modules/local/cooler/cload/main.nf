process COOLER_CLOAD {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::cooler=0.8.11"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'quay.io/biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(pairs), path(index)
    val cool_bin
    path chromsizes

    output:
    tuple val(meta), val(cool_bin), path("*.cool"), emit: cool
    path "versions.yml"                           , emit: versions

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args     = task.ext.args.tokenize()
    def tool     = (task.ext.args.contains('hiclib')) ? "hiclib" :
        (task.ext.args.contains('tabix')) ? 'tabix' :
        (task.ext.args.contains('pairs')) ? 'pairs' : 'pairix'
    def nproc    = tool in ['pairix','tabix'] ? "--nproc ${task.cpus}" : ''
    args.removeIf { it.contains('hiclib') }
    args.removeIf { it.contains('tabix') }
    args.removeIf { it.contains('pairs') }
    args.removeIf { it.contains('pairix') }

    """
    cooler cload $tool \\
        ${args.join(' ')} \\
        $nproc \\
        ${chromsizes}:${cool_bin} \\
        $pairs \\
        ${prefix}.${cool_bin}.cool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(echo \$(cooler --version 2>&1) | sed 's/cooler, version //')
    END_VERSIONS
    """
}
