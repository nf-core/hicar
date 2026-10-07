process HOMER_INSTALL {
    tag "$genome"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::homer=4.11 bioconda::samtools=1.11 conda-forge::r-base=4.0.2 bioconda::bioconductor-deseq2=1.30.0 bioconda::bioconductor-edger=3.32.0 anaconda::perl=5.26.2 wget"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0' :
        'biocontainers/mulled-v2-29293b111ffe5b4c1d1e14c711264aaed6b97b4a:594338b771cacf1623bd27772b5e12825f8835f2-0' }"

    input:
    val genome

    output:
    stdout  emit: output
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ## install genome
    homerDir=\$(grep "my \\\$homeDir =" \$(which runHiCpca.pl))
    homerDir=\${homerDir/my \\\$homeDir = /}
    homerDir=\${homerDir/;/}
    homerDir=\${homerDir//[\\"\\']}
    perl \${homerDir}configureHomer.pl -install $genome
    echo 'done'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: \$(echo \$(perl \${homerDir}configureHomer.pl -list 2>&1) | grep Executables | sed 's/^.*homer\\s\\+//; s/ Code.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
