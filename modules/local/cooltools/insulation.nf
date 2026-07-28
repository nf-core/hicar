process COOLTOOLS_INSULATION {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::cooltools=0.5.1 bioconda::ucsc-bedgraphtobigwig=377"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' :
        'biocontainers/mulled-v2-c81d8d6b6acf4714ffaae1a274527a41958443f6:cc7ea58b8cefc76bed985dcfe261cb276ed9e0cf-0' }"

    input:
    tuple val(meta), path(cool)
    val resolution

    output:
    tuple val(meta), path("*tsv")                , emit:results
    tuple val(meta), path("*.bed"), optional:true, emit:tads
    path("versions.yml")                         , emit:versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def res = [ 3*resolution,
                5*resolution,
                10*resolution,
                25*resolution].join(' ').trim()
    """
    nproc=\$(python -c "import cooltools;import sys;from packaging import version; print('-p '+sys.argv[1]) if version.parse(cooltools.__version__) >= version.parse('0.5.4') else print('')" $task.cpus)
    cooltools insulation \\
        \$nproc \\
        $args \\
        -o ${prefix}_insulation.tsv \\
        $cool \\
        $res

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | sed 's/cooltools, version //')
    END_VERSIONS
    """
}
