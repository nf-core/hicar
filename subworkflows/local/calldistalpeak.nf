/*
 * call peak by MACS2 for ATAC reads
 */
params.options = [:]

include { R1READS             } from '../../modules/local/atacreads/r1reads'          addParams(options: params.options.r1reads)
include { MERGEREADS
    as MERGE_R1READS          } from '../../modules/local/atacreads/mergereads'       addParams(options: params.options.merge_r1reads)
include { MACS2_CALLPEAK
    as MACS2_CALLR1PEAK       } from '../../modules/local/atacreads/macs2'            addParams(options: params.options.macs2_callr1peak)
include { DUMPREADS
    as DUMPR1READS            } from '../../modules/local/atacreads/dumpreads'        addParams(options: params.options.dump_r1_reads_per_group)
include { DUMPREADS
    as DUMPR1READS_SAMPLE     } from '../../modules/local/atacreads/dumpreads'        addParams(options: params.options.dump_r1_reads_per_sample)
include { MERGE_PEAK
    as MERGE_R1PEAK           } from '../../modules/local/atacreads/mergepeak'        addParams(options: params.options.merge_r1peak)
include { BEDTOOLS_GENOMECOV
    as BEDTOOLS_GENOMECOV_R1SAM } from '../../modules/nf-core/modules/bedtools/genomecov/main'  addParams(options: params.options.bedtools_genomecov_per_sample)
include { BEDTOOLS_SORT       } from '../../modules/nf-core/modules/bedtools/sort/main'  addParams(options: params.options.bedtools_sort_per_group)
include { BEDTOOLS_SORT
    as BEDTOOLS_SORT_R1SAM    } from '../../modules/nf-core/modules/bedtools/sort/main'  addParams(options: params.options.bedtools_sort_per_sample)
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
    ch_version = R1READS(distalpair).version

    // merge the read in same group
    R1READS.out.bed
            .map{meta, bed -> [meta.group, bed]}
            .groupTuple()
            .map{it -> [[id:it[0]], it[1]]} // id is group
            .set{read4merge}
    MERGE_R1READS(read4merge)
    ch_version = ch_version.mix(MERGE_R1READS.out.version)

    // call ATAC narrow peaks for group
    MACS2_CALLR1PEAK(MERGE_R1READS.out.bed, macs_gsize)
    ch_version = ch_version.mix(MACS2_CALLR1PEAK.out.version)

    // merge peaks
    r1_peaks = MACS2_CALLR1PEAK.out.peak.map{it[1]}.collect()
    MERGE_R1PEAK(r1_peaks)

    // stats
    //R1QC(atac_peaks, MERGE_R1READS.out.bed.map{it[1]}.collect(), gtf)
    //ch_version = ch_version.mix(R1QC.out.version)

    // dump R1 reads for each group for maps
    DUMPR1READS(MERGE_R1READS.out.bed)
    BEDTOOLS_SORT(MACS2_CALLR1PEAK.out.pileup)
    UCSC_BEDCLIP(BEDTOOLS_SORT.out.bed, chromsizes)
    UCSC_R1BEDGRAPHTOBIGWIG(UCSC_BEDCLIP.out.bedgraph, chromsizes)

    // dump ATAC reads for each samples for differential analysis
    DUMPR1READS_SAMPLE(R1READS.out.bed)
    ch_version = ch_version.mix(DUMPR1READS_SAMPLE.out.version)
    BEDTOOLS_GENOMECOV_R1SAM(R1READS.out.bed, chromsizes, "bedgraph")
    BEDTOOLS_SORT_R1SAM(BEDTOOLS_GENOMECOV_R1SAM.out.genomecov)
    ch_version = ch_version.mix(BEDTOOLS_GENOMECOV_R1SAM.out.version)
    UCSC_BEDGRAPHTOBIGWIG_R1SAM(BEDTOOLS_SORT_R1SAM.out.bed, chromsizes)
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_R1SAM.out.version)

    emit:
    peak       = MACS2_CALLR1PEAK.out.peak              // channel: [ val(meta), path(peak) ]
    xls        = MACS2_CALLR1PEAK.out.xls               // channel: [ val(meta), path(xls) ]
    mergedpeak = MERGE_R1PEAK.out.peak                  // channel: [ path(bed) ]
    //stats      = R1QC.out.stats                       // channel: [ path(csv) ]
    reads      = DUMPR1READS.out.peak                   // channel: [ val(meta), path(bedgraph) ]
    samplereads= DUMPR1READS.out.peak                   // channel: [ val(meta), path(bedgraph) ]
    bws        = UCSC_R1BEDGRAPHTOBIGWIG.out.bigwig     // channel: [ val(meta), path(bigwig) ]
    version    = ch_version                             // channel: [ path(version) ]
}