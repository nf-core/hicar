process PEAKACHU_MODEL {
    label 'process_single'

    conda "bioconda::cooler=0.8.11"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooler:0.8.11--pyh3252c3a_0' :
        'biocontainers/cooler:0.8.11--pyh3252c3a_0' }"

    input:
    tuple val(meta), path(cool)

    output:
    tuple val(meta), stdout                                , emit: model
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bin=\$((${meta.bin}/1000))
    reads=\$(cooler info -f sum $cool)
    reads=\$((\$reads/1000000))
    arr=(10 20 30 40 50 60 70 80 90 100 125 150 175 200 225 250 275 300)
    min=10000
    roundnum=10
    for i in \${arr[@]}
    do
        cur=\$((\$reads-\$i))
        cur=\${cur#-}
        if [ \$cur -lt \$min ]
        then
            min=\$cur
            roundnum=\$i
        else
            break
        fi
    done
    echo "${params.peakachu_pretrained_url}\${roundnum}million.\${bin}kb.pkl"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooler: \$(echo \$(cooler --version 2>&1) | sed 's/cooler, version //')
    END_VERSIONS
    """
}
