/*
 * Uncompress and prepare reference genome files
 */

include {
    GUNZIP as GUNZIP_FASTA;
    GUNZIP as GUNZIP_GTF;
    GUNZIP as GUNZIP_GFF;
    GUNZIP as GUNZIP_GENE_BED;
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/gunzip/main'
include { GTF2BED                     } from '../../modules/local/gtf2bed'
include { CHROMSIZES                  } from '../../modules/local/genome/chromsizes'
include { GENOME_FILTER               } from '../../modules/local/genome/filter'
include { COOLER_DIGEST               } from '../../modules/nf-core/cooler/digest/main'
include { RE_CUTSITE                  } from '../../modules/local/re_cut'
include { GFFREAD                     } from '../../modules/nf-core/gffread/main'
include { GENMAP_INDEX                } from '../../modules/nf-core/genmap/index/main'
include { GENMAP_MAPPABILITY          } from '../../modules/nf-core/genmap/mappability/main'
include { UCSC_WIGTOBIGWIG            } from '../../modules/nf-core/ucsc/wigtobigwig/main'
include { BWA_INDEX                   } from '../../modules/nf-core/bwa/index/main'

workflow PREPARE_GENOME {
    main:

    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA ( [[id:'fasta'], file("${params.fasta}", checkIfExists: true)] )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map{it[1]}
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
            GUNZIP_GTF ( [[id:'gtf'], file("${params.gtf}", checkIfExists: true)] )
            ch_gtf = GUNZIP_GTF.out.gunzip.map{it[1]}
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF ( [[id:'gff'], file("${params.gff}", checkIfExists: true)] )
            ch_gff = GUNZIP_GFF.out.gunzip.map{it[1]}
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_version = ch_version.mix(GFFREAD.out.versions)
    }

    /*
     * Uncompress gene BED annotation file or create from GTF if required
     */
    ch_gene_bed = Channel.empty()
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            GUNZIP_GENE_BED ( [[id:'gene_bed'], file("${params.gene_bed}", checkIfExists: true)] )
            ch_gene_bed = GUNZIP_GENE_BED.out.gunzip.map{it[1]}
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_version = ch_version.mix(GTF2BED.out.versions)
    }

    /*
     * Create chromosome sizes file
     */
    ch_chrom_sizes = CHROMSIZES ( ch_fasta ).sizes

    /*
     * Calculate effective genome sizes
     */
    genome_size = 0
    if (params.macs_gsize) {
        genome_size = params.macs_gsize
    } else {
        //genome size remove all N
        // ref:https://github.com/macs3-project/MACS/issues/299
        // https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html method 1
        genome_size = ch_fasta.splitFasta( record: [id: true, seqString: true] )
            .map { it.seqString.replaceAll('[^acgtACGT]','').length() }
            .sum()
    }

    /*
     * Prepare ucsc annotation name
     */
    ucscname = null
    if (params.ucscname) {
        ucscname = params.ucscname
    } else {
        ucsc_map = ["GRCh38":"hg38", "GRCh37":"hg19",
                    "GRCm38":"mm10", "TAIR10":"tair10",
                    "UMD3.1":"bosTau8", "CanFam3.1":"canFam3",
                    "WBcel235":"ce11", "GRCz10":"danRer10",
                    "BDGP6":"dm6", "EquCab2":'equCab2',
                    "Galgal4":"galGal4", "CHIMP2.1.4":"panTro4",
                    "Rnor_6.0":"rn6", "R64-1-1":"sacCer3",
                    "Sscrofa10.2":"susScr3"]
        if(params.genome){
            if(ucsc_map[params.genome]){
                ucscname = ucsc_map[params.genome]
            }else{
                ucscname = params.genome
            }
        }
    }

    if (params.blacklist) {
        ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true)
    } else { ch_blacklist = Channel.empty() }
    ch_blacklist = ch_blacklist.ifEmpty([])

    filtered_bed = GENOME_FILTER (
        ch_chrom_sizes,
        ch_blacklist
    ).bed
    ch_version = ch_version.mix(GENOME_FILTER.out.versions)

    /*
     * Create digest genome file for PAIRTOOLS_PAIRE
     */
    digest_genome_bed = COOLER_DIGEST (
        ch_fasta,
        ch_chrom_sizes,
        params.enzyme
    ).bed
    ch_version = ch_version.mix(COOLER_DIGEST.out.versions)

    /*
     * get enzyme cut site and position for function maps:cut or enzyme_cut
     */
    RE_CUTSITE ( params.enzyme )
    ch_version = ch_version.mix(RE_CUTSITE.out.versions)

    /*
     * mappability
     */
    if(!params.mappability){
        GENMAP_INDEX(ch_fasta).index | GENMAP_MAPPABILITY
        ch_version = ch_version.mix(GENMAP_MAPPABILITY.out.versions)
        ch_mappability = UCSC_WIGTOBIGWIG(
            GENMAP_MAPPABILITY.out.wig.map{[[id:'mappability'], it]},
            ch_chrom_sizes).bw.map{it[1]}
        ch_version = ch_version.mix(UCSC_WIGTOBIGWIG.out.versions)
    }else{
        ch_mappability = Channel.fromPath(params.mappability, checkIfExists: true)
    }

    /*
     * Uncompress bwa index or generate from scratch if required
     */
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)).map{ [['id':ucscname], it] } : BWA_INDEX ( ch_fasta.map{[['id':ucscname], it]} ).index

    emit:
    fasta             = ch_fasta                       // path: genome.fasta,
    gtf               = ch_gtf                         // path: genome.gtf,
    gene_bed          = ch_gene_bed                    // path: gene.bed,
    chrom_sizes       = ch_chrom_sizes                 // path: genome.sizes,
    blacklist         = ch_blacklist                   // path: blacklist.bed,
    bed               = filtered_bed                   // path: *.bed,
    digest_genome     = digest_genome_bed              // path: bed
    site              = RE_CUTSITE.out.site            // value: site 5position
    mappability       = ch_mappability                 // path: bw
    bwa_index         = ch_bwa_index                   // path: bwt,amb,sa,ann,pac
    gsize             = genome_size                    // value: macs2 genome size
    ucscname          = ucscname                       // value: ucsc annotation name
    versions          = ch_version                     // path: *.version.yml
}
