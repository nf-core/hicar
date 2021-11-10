// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process READS_STAT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "r::r-magrittr=1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/r-magrittr:1.5--r3.2.2_0"
    } else {
        container "quay.io/biocontainers/r-magrittr:1.5--r3.2.2_0"
    }

    input:
    tuple val(meta), path(raw), path(dedup)

    output:
    tuple val(meta), path("*.csv"), emit: stat
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    #!/usr/bin/env Rscript
    ## generate the statistis for each samples
    versions <- c("${getProcessName(task.process)}:",
        paste0("    R:", paste(R.version\$major, R.version\$minor, sep=".")))
    writeLines(versions, "versions.yml") # wirte versions.yml

    raw = "$raw"
    dedup = "$dedup"
    out = "${prefix}.reads_stats.csv"

    sample_name <- sub(".raw.pairsam.stat", "", basename(raw))
    getDat <- function(f){
        dat <- read.delim(f, header=FALSE)
        res <- dat[, 2]
        names(res) <- dat[, 1]
        return(res)
    }
    all_pairs <- getDat(raw)
    dep_pairs <- getDat(dedup)

    df <- data.frame(sample=sample_name,
            total=all_pairs["total"],
            duplicate=dep_pairs['total_dups'],
            non_duplicated=dep_pairs['total_nodups'],
            duplication_rate=round(100*dep_pairs['total_dups']/all_pairs["total"],2),
            trans=dep_pairs['trans'],
            cis=dep_pairs['cis'],
            longRange=dep_pairs['cis_20kb+'])
    write.csv(df, out, row.names=FALSE)
    """
}
