/*
 * call peak by MACS2 for ATAC reads
 */
include { PAIRTOOLS_SELECT
    as  PAIRTOOLS_SELECT_SHORT} from '../../modules/nf-core/modules/pairtools/select/main'
include { SHIFT_READS         } from '../../modules/local/atacreads/shiftreads'
include { MERGE_READS         } from '../../modules/local/atacreads/mergereads'
include { MACS2_CALLPEAK      } from '../../modules/nf-core/modules/macs2/callpeak/main'
include { DUMP_READS          } from '../../modules/local/atacreads/dumpreads'
include { DUMP_READS
    as DUMP_READS_PER_SAMPLE  } from '../../modules/local/atacreads/dumpreads'
include { MERGE_PEAK          } from '../../modules/local/atacreads/mergepeak'
include { ATACQC              } from '../../modules/local/atacreads/atacqc'
include { BEDTOOLS_GENOMECOV
    as BEDTOOLS_GENOMECOV_PER_SAMPLE } from '../../modules/nf-core/modules/bedtools/genomecov/main'
include { BEDFILES_SORT
    as BEDFILES_SORT_PER_GROUP       } from '../../modules/local/atacreads/bedsort'
include { BEDFILES_SORT
    as BEDFILES_SORT_PER_SAMPLE      } from '../../modules/local/atacreads/bedsort'
include { UCSC_BEDCLIP        } from '../../modules/nf-core/modules/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_PER_GROUP} from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE} from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'

workflow ATAC_PEAK {
    take:
    validpair  // channel: [ val(meta), [pairs] ]
    chromsizes // channel: [ path(size) ]
    macs_gsize // channel: value
    gtf        // channel: [ path(gtf) ]

    main:
    // extract ATAC reads, split the pairs into longRange_Trans pairs and short pairs
    ch_version = PAIRTOOLS_SELECT_SHORT(validpair).versions
    // shift Tn5 insertion for longRange_Trans pairs
    SHIFT_READS(PAIRTOOLS_SELECT_SHORT.out.unselected)
    ch_version = ch_version.mix(SHIFT_READS.out.versions)

    // merge the read in same group
    SHIFT_READS.out.bed
            .map{meta, bed -> [meta.group, bed]}
            .groupTuple()
            .map{it -> [[id:it[0]], it[1]]} // id is group
            .set{read4merge}
    MERGE_READS(read4merge)
    ch_version = ch_version.mix(MERGE_READS.out.versions)

    // call ATAC narrow peaks for group
    MACS2_CALLPEAK(MERGE_READS.out.bed.map{[it[0], it[1], []]}, macs_gsize)
    ch_version = ch_version.mix(MACS2_CALLPEAK.out.versions)

    // merge peaks
    atac_peaks = MACS2_CALLPEAK.out.peak.map{it[1]}.collect()
    MERGE_PEAK(atac_peaks)

    // stats
    ATACQC(atac_peaks, MERGE_READS.out.bed.map{it[1]}.collect(), gtf)
    ch_version = ch_version.mix(ATACQC.out.versions)

    // dump ATAC reads for each group for maps
    DUMP_READS(MERGE_READS.out.bed)
    BEDFILES_SORT_PER_GROUP(MACS2_CALLPEAK.out.bdg.map{[it[0], it[1].findAll{it.toString().contains('pileup')}]}, "bedgraph")
    UCSC_BEDCLIP(BEDFILES_SORT_PER_GROUP.out.sorted, chromsizes)
    UCSC_BEDGRAPHTOBIGWIG_PER_GROUP(UCSC_BEDCLIP.out.bedgraph, chromsizes)

    // dump ATAC reads for each samples for differential analysis
    DUMP_READS_PER_SAMPLE(SHIFT_READS.out.bed)
    ch_version = ch_version.mix(DUMP_READS.out.versions)
    BEDTOOLS_GENOMECOV_PER_SAMPLE(SHIFT_READS.out.bed.map{[it[0], it[1], "1"]}, chromsizes, "bedgraph")
    BEDFILES_SORT_PER_SAMPLE(BEDTOOLS_GENOMECOV_PER_SAMPLE.out.genomecov, "bedgraph")
    ch_version = ch_version.mix(BEDTOOLS_GENOMECOV_PER_SAMPLE.out.versions)
    UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE(BEDFILES_SORT_PER_SAMPLE.out.sorted, chromsizes)
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_PER_SAMPLE.out.versions)

    emit:
    peak       = MACS2_CALLPEAK.out.peak              // channel: [ val(meta), path(peak) ]
    xls        = MACS2_CALLPEAK.out.xls               // channel: [ val(meta), path(xls) ]
    mergedpeak = MERGE_PEAK.out.peak                  // channel: [ path(bed) ]
    stats      = ATACQC.out.stats                     // channel: [ path(csv) ]
    reads      = DUMP_READS.out.peak                  // channel: [ val(meta), path(bedgraph) ]
    samplereads= DUMP_READS.out.peak                  // channel: [ val(meta), path(bedgraph) ]
    bws        = UCSC_BEDGRAPHTOBIGWIG_PER_GROUP.out.bigwig     // channel: [ val(meta), path(bigwig) ]
    versions   = ch_version                           // channel: [ path(version) ]
}
