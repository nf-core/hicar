/*
 * call peak by MACS2 for ATAC reads
 */
params.options = [:]

include { R1READS             } from '../../modules/local/atacreads/r1reads'          addParams(options: params.options.r1reads)
include { MERGEREADS
    as MERGE_R1READS          } from '../../modules/local/atacreads/mergereads'       addParams(options: params.options.merge_r1reads)
include { MACS2_CALLPEAK
    as MACS2_CALLR1PEAK       } from '../../modules/nf-core/modules/macs2/callpeak/main'            addParams(options: params.options.macs2_callr1peak)
include { DUMPREADS
    as DUMPR1READS            } from '../../modules/local/atacreads/dumpreads'        addParams(options: params.options.dump_r1_reads_per_group)
include { DUMPREADS
    as DUMPR1READS_SAMPLE     } from '../../modules/local/atacreads/dumpreads'        addParams(options: params.options.dump_r1_reads_per_sample)
include { MERGE_PEAK
    as MERGE_R1PEAK           } from '../../modules/local/atacreads/mergepeak'        addParams(options: params.options.merge_r1peak)
include { BEDTOOLS_GENOMECOV
    as BEDTOOLS_GENOMECOV_R1SAM } from '../../modules/nf-core/modules/bedtools/genomecov/main'  addParams(options: params.options.bedtools_genomecov_per_sample)
include { BEDFILES_SORT       } from '../../modules/local/atacreads/bedsort'          addParams(options: params.options.bedtools_sort_per_group)
include { BEDFILES_SORT
    as BEDFILES_SORT_SAM      } from '../../modules/local/atacreads/bedsort'          addParams(options: params.options.bedtools_sort_per_sample)
include { UCSC_BEDCLIP        } from '../../modules/nf-core/modules/ucsc/bedclip/main'  addParams(options: params.options.ucsc_bedclip)
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_R1BEDGRAPHTOBIGWIG     } from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'  addParams(options: params.options.ucsc_bedgraphtobigwig_per_r1_group)
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_R1SAM } from '../../modules/nf-core/modules/ucsc/bedgraphtobigwig/main'  addParams(options: params.options.ucsc_bedgraphtobigwig_per_r1_sample)

workflow R1_PEAK {
    take:
    distalpair  // channel: [ val(meta), [pairs] ]
    chromsizes // channel: [ path(size) ]
    macs_gsize // channel: value
    gtf        // channel: [ path(gtf) ]

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

    // call ATAC narrow peaks for group
    MACS2_CALLR1PEAK(MERGE_R1READS.out.bed.map{[it[0], it[1], []]}, macs_gsize)
    ch_version = ch_version.mix(MACS2_CALLR1PEAK.out.versions)

    // merge peaks
    r1_peaks = MACS2_CALLR1PEAK.out.peak.map{it[1]}.collect()
    MERGE_R1PEAK(r1_peaks)

    // stats
    //R1QC(atac_peaks, MERGE_R1READS.out.bed.map{it[1]}.collect(), gtf)
    //ch_version = ch_version.mix(R1QC.out.versions)

    // dump R1 reads for each group for maps
    DUMPR1READS(MERGE_R1READS.out.bed)
    BEDFILES_SORT(MACS2_CALLR1PEAK.out.bdg.map{[it[0], it[1].findAll{it.toString().contains('pileup')}]}, "bedgraph")
    UCSC_BEDCLIP(BEDFILES_SORT.out.sorted, chromsizes)
    UCSC_R1BEDGRAPHTOBIGWIG(UCSC_BEDCLIP.out.bedgraph, chromsizes)

    // dump ATAC reads for each samples for differential analysis
    DUMPR1READS_SAMPLE(R1READS.out.bed)
    ch_version = ch_version.mix(DUMPR1READS_SAMPLE.out.versions)
    BEDTOOLS_GENOMECOV_R1SAM(R1READS.out.bed.map{[it[0], it[1], "1"]}, chromsizes, "bedgraph")
    BEDFILES_SORT_SAM(BEDTOOLS_GENOMECOV_R1SAM.out.genomecov, "bedgraph")
    ch_version = ch_version.mix(BEDTOOLS_GENOMECOV_R1SAM.out.versions)
    UCSC_BEDGRAPHTOBIGWIG_R1SAM(BEDFILES_SORT_SAM.out.sorted, chromsizes)
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_R1SAM.out.versions)

    emit:
    peak       = MACS2_CALLR1PEAK.out.peak              // channel: [ val(meta), path(peak) ]
    xls        = MACS2_CALLR1PEAK.out.xls               // channel: [ val(meta), path(xls) ]
    mergedpeak = MERGE_R1PEAK.out.peak                  // channel: [ path(bed) ]
    //stats      = R1QC.out.stats                       // channel: [ path(csv) ]
    reads      = DUMPR1READS.out.peak                   // channel: [ val(meta), path(bedgraph) ]
    samplereads= DUMPR1READS.out.peak                   // channel: [ val(meta), path(bedgraph) ]
    bws        = UCSC_R1BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), path(bigwig) ]
    versions   = ch_version                             // channel: [ path(version) ]
}
