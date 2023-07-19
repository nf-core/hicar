process GENOME_FILTER {
    label 'process_low'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    path sizes
    path blacklist

    output:
    path "*.bed",         emit: bed
    path "versions.yml",  emit: versions

    script:
    def file_out = "${sizes.simpleName}.include_regions.bed"
    if (blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bedtools: \$(echo \$(bedtools --version) | sed -e "s/bedtools v//g")
        END_VERSIONS
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            awk: \$(echo \$(awk --version 2>&1 || awk -W version 2>&1) | sed 's/[[:alpha:]|(|)|[:space:]]//g; s/,.*\$//')
        END_VERSIONS
        """
    }
}
