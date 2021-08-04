/*
 * Uncompress and prepare reference genome files
 */
params.options = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/modules/gunzip/main'               addParams( options: params.options.gunzip        )
include { GTF2BED                     } from '../../modules/local/gtf2bed'                             addParams( options: params.options.gtf2bed       )
include { CHROMSIZES                  } from '../../modules/local/genome/chromsizes'                   addParams( options: params.options.chromsizes    )
include { GENOME_FILTER               } from '../../modules/local/genome/filter'                       addParams( options: params.options.genomefilter  )
include { COOLER_DIGEST               } from '../../modules/nf-core/modules/cooler/digest/main'        addParams( options: params.options.digest_genome )
include { GFFREAD                     } from '../../modules/nf-core/modules/gffread/main'              addParams( options: params.options.gffread       )
include { BWA_INDEX                   } from '../../modules/nf-core/modules/bwa/index/main'            addParams( options: params.options.bwa_index     )

workflow PREPARE_GENOME {
    main:

    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * Uncompress GTF annotation file or create from GFF3 if required
     */
    ch_version = Channel.empty()
    ch_gtf = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_version = ch_version.mix(GFFREAD.out.version.ifEmpty(null))
    }

    /*
     * Uncompress gene BED annotation file or create from GTF if required
     */
    ch_gene_bed = Channel.empty()
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( params.gene_bed ).gunzip
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_version = ch_version.mix(GTF2BED.out.version.ifEmpty(null))
    }

    /*
     * Create chromosome sizes file
     */
    ch_chrom_sizes = CHROMSIZES ( ch_fasta ).sizes

    if (params.blacklist) {
        ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true)
    } else { ch_blacklist = Channel.empty() }
    ch_blacklist = ch_blacklist.ifEmpty([])

    filtered_bed = GENOME_FILTER (
        ch_chrom_sizes,
        ch_blacklist
    ).bed
    ch_version = ch_version.mix(GENOME_FILTER.out.version.ifEmpty(null))

    /*
     * Create digest genome file for PAIRTOOLS_PAIRE
     */
     digest_genome_bed = COOLER_DIGEST (
        ch_fasta,
        ch_chrom_sizes,
        params.enzyme
    ).bed
    ch_version = ch_version.mix(COOLER_DIGEST.out.version.ifEmpty(null))

    /*
     * Uncompress bwa index or generate from scratch if required
     */
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta ).index

    emit:
    fasta             = ch_fasta                       // path: genome.fasta,
    gtf               = ch_gtf                         // path: genome.gtf,
    gene_bed          = ch_gene_bed                    //           path: gene.bed,
    chrom_sizes       = ch_chrom_sizes                 //           path: genome.sizes,
    blacklist         = ch_blacklist                   //           path: blacklist.bed,
    bed               = filtered_bed                   //           path: *.bed,
    digest_genome     = digest_genome_bed              //           path: bed
    bwa_index         = ch_bwa_index                   //           path: fasta
    version           = ch_version                       // path: *.version.txt
}
