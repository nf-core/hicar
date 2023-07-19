/*
 * call peak by fragment reads
 */
include { R1READS             } from '../../../modules/local/fragmentreads/r1reads'
include { MERGE_READS
    as MERGE_R1READS          } from '../../../modules/local/atacreads/mergereads'
include { CALL_R1PEAK         } from '../../../modules/local/fragmentreads/call_peak'
include {
    DUMP_READS
        as DUMP_R1_READS_PER_GROUP;
    DUMP_READS
        as DUMP_R1_READS_PER_SAMPLE    } from '../../../modules/local/atacreads/dumpreads'
include { MERGE_PEAK
    as MERGE_R1PEAK           } from '../../../modules/local/atacreads/mergepeak'
include { BEDTOOLS_GENOMECOV
    as BEDTOOLS_GENOMECOV_PER_R1_SAMPLE } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDFILES_SORT
    as BEDFILES_SORT_PER_GROUP       } from '../../../modules/local/atacreads/bedsort'
include { BEDFILES_SORT
    as BEDFILES_SORT_PER_SAMPLE      } from '../../../modules/local/atacreads/bedsort'
include { UCSC_BEDCLIP        } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_PER_R1_GROUP  } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_PER_R1_SAMPLE } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow R1_PEAK {
    take:
    distalpair // channel: [ val(meta), [pairs] ]
    chromsizes // channel: [ path(size) ]
    cut        // channel: [ path(cut) ]
    gtf        // channel: [ path(gtf) ]
    pval       // val
    short_bed_postfix

    main:
    // extract and sort R1 reads
    ch_version = R1READS(distalpair).versions

    // merge the read in same group
    R1READS.out.bed
            .map{meta, bed -> [meta.group, bed]}
            .groupTuple()
            .map{it -> [[id:it[0]], it[1]]} // id is group
            .set{read4merge}
    MERGE_R1READS(read4merge)
    ch_version = ch_version.mix(MERGE_R1READS.out.versions)

    // call fragment narrow peaks for group
    CALL_R1PEAK(MERGE_R1READS.out.bed, cut, pval)
    ch_version = ch_version.mix(CALL_R1PEAK.out.versions)

    // merge peaks
    r1_peaks = CALL_R1PEAK.out.peak.map{it[1]}.collect()
    MERGE_R1PEAK(r1_peaks)

    // dump R1 reads for each group for maps
    DUMP_R1_READS_PER_GROUP(MERGE_R1READS.out.bed, short_bed_postfix)
    BEDFILES_SORT_PER_GROUP(CALL_R1PEAK.out.bdg, "bedgraph")
    UCSC_BEDCLIP(BEDFILES_SORT_PER_GROUP.out.sorted, chromsizes)
    UCSC_BEDGRAPHTOBIGWIG_PER_R1_GROUP(UCSC_BEDCLIP.out.bedgraph, chromsizes)

    // dump R1 reads for each samples for differential analysis
    DUMP_R1_READS_PER_SAMPLE(R1READS.out.bed, short_bed_postfix)
    ch_version = ch_version.mix(DUMP_R1_READS_PER_SAMPLE.out.versions)
    BEDTOOLS_GENOMECOV_PER_R1_SAMPLE(R1READS.out.bed.map{[it[0], it[1], "1"]}, chromsizes, "bedgraph")
    BEDFILES_SORT_PER_SAMPLE(BEDTOOLS_GENOMECOV_PER_R1_SAMPLE.out.genomecov, "bedgraph")
    ch_version = ch_version.mix(BEDTOOLS_GENOMECOV_PER_R1_SAMPLE.out.versions)
    UCSC_BEDGRAPHTOBIGWIG_PER_R1_SAMPLE(BEDFILES_SORT_PER_SAMPLE.out.sorted, chromsizes)
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_PER_R1_SAMPLE.out.versions)

    emit:
    peak       = CALL_R1PEAK.out.peak                   // channel: [ val(meta), path(peak) ]
    mergedpeak = MERGE_R1PEAK.out.peak                  // channel: [ path(bed) ]
    //stats      = R1QC.out.stats                       // channel: [ path(csv) ]
    reads      = DUMP_R1_READS_PER_GROUP.out.peak       // channel: [ val(meta), path(bedgraph) ]
    samplereads= DUMP_R1_READS_PER_GROUP.out.peak       // channel: [ val(meta), path(bedgraph) ]
    bws        = UCSC_BEDGRAPHTOBIGWIG_PER_R1_GROUP.out.bigwig     // channel: [ val(meta), path(bigwig) ]
    versions   = ch_version                             // channel: [ path(version) ]
}
