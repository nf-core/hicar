process VIRTUAL4C_BY_COOLTOOLS {
    tag "${bin_size}"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::cooltools=0.5.4"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cooltools:0.5.4--py39hbf8eff0_0' :
        'biocontainers/cooltools:0.5.4--py39hbf8eff0_0' }"

    input:
    tuple val(bin_size), path(cool), path(viewpoint)
    val v4c_max_events

    output:
    path "${prefix}/*"            , emit: v4c
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${bin_size}"
    """
    # VIEWPOINT : the viewpoint to use for the virtual 4C profile. Provide as a UCSC-string (e.g. chr1:1-1000)
    mkdir -p $prefix
    for vpfile in $viewpoint
    do
        evt=0
        cat \$vpfile | while read -r vp
        do
            if [ \$evt -gt $v4c_max_events ]; then
                break
            fi
            vp=(\${vp//[^a-zA-Z0-9]/ })
            if [[ \${vp[1]} =~ ^[0-9]+\$ ]]; then
                evt=\$((\$evt+1))
                vp=\${vp[0]}:\${vp[1]}-\${vp[2]}
                for cool_file in $cool
                do
                    cooltools virtual4c \\
                        $args \\
                        -p ${task.cpus} \\
                        -o ${prefix}/\${cool_file}_\$vp \\
                        \$cool_file \\
                        \$vp
                done
            fi
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cooltools: \$(cooltools --version 2>&1 | sed 's/cooletools, version //')
    END_VERSIONS
    """
}
