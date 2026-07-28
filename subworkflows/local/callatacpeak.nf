/*
 * call peak by MACS2 for ATAC reads
 */
include { PAIRTOOLS_SELECT
    as  PAIRTOOLS_SELECT_SHORT              } from '../../modules/nf-core/pairtools/select/main'
include { SHIFT_READS                       } from '../../modules/local/atacreads/shiftreads'
include { MERGE_READS                       } from '../../modules/local/atacreads/mergereads'
include { MACS2_CALLPEAK                    } from '../../modules/nf-core/macs2/callpeak/main'
include {
    DUMP_READS
        as DUMP_READS_PER_GROUP;
    DUMP_READS
        as DUMP_READS_PER_SAMPLE            } from '../../modules/local/atacreads/dumpreads'
include { MERGE_PEAK                        } from '../../modules/local/atacreads/mergepeak'
include { ATACQC                            } from '../../modules/local/atacreads/atacqc'
include { COVERAGE_SCALE                    } from './coveragescale'
include {
    BEDTOOLS_GENOMECOV
        as BEDTOOLS_GENOMECOV_PER_SAMPLE;
    BEDTOOLS_GENOMECOV
        as BEDTOOLS_GENOMECOV_PER_GROUP     } from '../../modules/nf-core/bedtools/genomecov/main'
include {
    BEDFILES_SORT
        as BEDFILES_SORT_PER_GROUP;
    BEDFILES_SORT
        as BEDFILES_SORT_PER_SAMPLE         } from '../../modules/local/atacreads/bedsort'
include { UCSC_BEDCLIP                      } from '../../modules/nf-core/ucsc/bedclip/main'
include {
    UCSC_BEDGRAPHTOBIGWIG
        as UCSC_BEDGRAPHTOBIGWIG_PER_GROUP;
    UCSC_BEDGRAPHTOBIGWIG
        as UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow ATAC_PEAK {
    take:
    validpair  // channel: [ val(meta), [pairs] ]
    chromsizes // channel: [ path(size) ]
    macs_gsize // channel: value
    gtf        // channel: [ path(gtf) ]
    method     // channel: value
    peak_file  // channel: value
    anchors    // channel: [ path(anchor_peaks) ]
    short_bed_postfix

    main:
    // extract ATAC reads, split the pairs into longRange_Trans pairs and short pairs
    ch_version = PAIRTOOLS_SELECT_SHORT(validpair).versions
    def doShift = method.toLowerCase() == "hicar"
    // shift Tn5 insertion for longRange_Trans pairs for hicar R2 reads
    SHIFT_READS(PAIRTOOLS_SELECT_SHORT.out.unselected, doShift)
    ch_version = ch_version.mix(SHIFT_READS.out.versions)

    // merge the read in same group
    SHIFT_READS.out.bed
            .map{meta, bed -> [meta.group, bed]}
            .groupTuple()
            .map{it -> [[id:it[0]], it[1]]} // id is group
            .set{read4merge}
    MERGE_READS(read4merge)
    ch_version = ch_version.mix(MERGE_READS.out.versions)

    if(peak_file){
        // use user defined peaks, and create bedgraph file for merged reads
        BEDTOOLS_GENOMECOV_PER_GROUP(MERGE_READS.out.bed.map{[it[0], it[1], "1"]}, chromsizes, "bedgraph")
        ch_group_bdg = BEDTOOLS_GENOMECOV_PER_GROUP.out.genomecov
        ch_group_peak = MERGE_READS.out.bed.map{[it[0], anchors]}
        anchor_peaks = anchors
    }else{
        // call ATAC narrow peaks for group
        MACS2_CALLPEAK(MERGE_READS.out.bed.map{[it[0], it[1], []]}, macs_gsize)
        ch_version = ch_version.mix(MACS2_CALLPEAK.out.versions)
        ch_group_bdg = MACS2_CALLPEAK.out.bdg.map{[it[0], it[1].findAll{it.toString().contains('pileup')}]}
        ch_group_peak= MACS2_CALLPEAK.out.peak
        anchor_peaks = MACS2_CALLPEAK.out.peak.map{it[1]}.collect()
    }

    // merge peaks
    MERGE_PEAK(anchor_peaks)

    // stats
    ATACQC(anchor_peaks, MERGE_READS.out.bed.map{it[1]}.collect(), gtf)
    ch_version = ch_version.mix(ATACQC.out.versions)

    // dump ATAC reads for each group for maps
    DUMP_READS_PER_GROUP(MERGE_READS.out.bed, short_bed_postfix)
    ch_version = ch_version.mix(DUMP_READS_PER_GROUP.out.versions)
    BEDFILES_SORT_PER_GROUP(ch_group_bdg, "bedgraph")
    UCSC_BEDCLIP(BEDFILES_SORT_PER_GROUP.out.sorted, chromsizes)
    UCSC_BEDGRAPHTOBIGWIG_PER_GROUP(UCSC_BEDCLIP.out.bedgraph, chromsizes)

    // dump ATAC reads for each samples for differential analysis
    DUMP_READS_PER_SAMPLE(SHIFT_READS.out.bed, short_bed_postfix)
    COVERAGE_SCALE(SHIFT_READS.out.counts)
    BEDTOOLS_GENOMECOV_PER_SAMPLE(
        SHIFT_READS.out.bed.join(COVERAGE_SCALE.out.scale),
        chromsizes,
        "bedgraph")
    BEDFILES_SORT_PER_SAMPLE(BEDTOOLS_GENOMECOV_PER_SAMPLE.out.genomecov, "bedgraph")
    ch_version = ch_version.mix(BEDTOOLS_GENOMECOV_PER_SAMPLE.out.versions)
    UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE(BEDFILES_SORT_PER_SAMPLE.out.sorted, chromsizes)
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE.out.versions)

    emit:
    peak       = ch_group_peak                        // channel: [ val(meta), path(peak) ]
    mergedpeak = MERGE_PEAK.out.peak                  // channel: [ path(bed) ]
    stats      = ATACQC.out.stats                     // channel: [ path(csv) ]
    reads      = DUMP_READS_PER_GROUP.out.peak        // channel: [ val(meta), path(bedgraph) ]
    samplereads= DUMP_READS_PER_SAMPLE.out.peak       // channel: [ val(meta), path(bedgraph) ]
    bws        = UCSC_BEDGRAPHTOBIGWIG_PER_GROUP.out.bigwig     // channel: [ val(meta), path(bigwig) ]
    versions   = ch_version                           // channel: [ path(version) ]
}
