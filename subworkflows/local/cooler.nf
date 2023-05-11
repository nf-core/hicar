/*
 * pair the proper mapped pairs
 */

include { COOLER_BALANCE } from '../../modules/nf-core/cooler/balance/main'
include { COOLER_CLOAD   } from '../../modules/nf-core/cooler/cload/main'
include { COOLER_MERGE   } from '../../modules/nf-core/cooler/merge/main'
include { COOLER_ZOOMIFY } from '../../modules/nf-core/cooler/zoomify/main'
include { COOLER_DUMP
    as COOLER_DUMP_PER_GROUP    } from '../../modules/nf-core/cooler/dump/main'
include { COOLER_DUMP
    as COOLER_DUMP_PER_SAMPLE   } from '../../modules/nf-core/cooler/dump/main'
include { DUMPREADS
    as DUMPREADS_PER_GROUP } from '../../modules/local/cooler/dumpreads'
include { DUMPREADS
    as DUMPREADS_PER_SAMPLE} from '../../modules/local/cooler/dumpreads'
include { JUICER_PRE       } from '../../modules/local/juicer/pre'

workflow COOLER {
    take:
    valid_pairs               // channel: [ val(meta), val(bin), [pairs], [pairs.px] ]
    chromsizes                // channel: [ path(chromsizes) ]
    hic_tools_jar             // channel: [ path(HiCTools jar for Pre) ]
    long_bedpe_postfix

    main:
    // HiC-like contact matrix
    ch_version = COOLER_CLOAD(valid_pairs.map{[it[0], it[2], it[3], it[1]]}, chromsizes).versions
    // Merge contacts
    COOLER_CLOAD.out.cool
                .map{
                    meta, cool, bin ->
                    [meta.group, bin, cool]
                }
                .groupTuple(by:[0, 1])
                .map{group, bin, cool -> [[id:group, bin:bin], cool]}
                .set{ch_cooler}
    COOLER_MERGE(ch_cooler)
    // create a balanced matrix for compartment and tad calls, see https://github.com/open2c/cooler/issues/48
    COOLER_BALANCE(COOLER_MERGE.out.cool.map{[it[0], it[1], false]})
    // create mcooler file for visualization
    COOLER_ZOOMIFY(COOLER_BALANCE.out.cool)
    // dump interaction bedpe for each group
    COOLER_DUMP_PER_GROUP(COOLER_MERGE.out.cool.map{[it[0], it[1], []]})
    // dump long interaction bedpe for each group for MAPS to call peaks
    DUMPREADS_PER_GROUP(COOLER_DUMP_PER_GROUP.out.bedpe, long_bedpe_postfix)
    ch_hic = Channel.empty()
    JUICER_PRE(DUMPREADS_PER_GROUP.out.gi, hic_tools_jar, chromsizes)
    ch_hic = JUICER_PRE.out.hic
    ch_version = ch_version.mix(JUICER_PRE.out.versions)

    // dump long interaction bedpe for each sample
    COOLER_DUMP_PER_SAMPLE(COOLER_CLOAD.out.cool.map{ meta, cool, bin -> [[id:meta.id, group:meta.group, bin:bin], cool, []]})
    DUMPREADS_PER_SAMPLE(COOLER_DUMP_PER_SAMPLE.out.bedpe, long_bedpe_postfix)
    ch_version = ch_version.mix(DUMPREADS_PER_GROUP.out.versions)

    emit:
    cool        = COOLER_BALANCE.out.cool                   // channel: [ val(meta), [cool] ]
    raw         = COOLER_MERGE.out.cool                     // channel: [ val(meta), [cool] ]
    mcool       = COOLER_ZOOMIFY.out.mcool                  // channel: [ val(meta), [mcool] ]
    hic         = ch_hic                                    // channel: [ val(meta), [hic] ]
    groupbedpe  = COOLER_DUMP_PER_GROUP.out.bedpe           // channel: [ val(meta), [bedpe] ]
    bedpe       = DUMPREADS_PER_GROUP.out.bedpe             // channel: [ val(meta), [bedpe] ]
    samplebedpe = DUMPREADS_PER_SAMPLE.out.bedpe            // channel: [ val(meta), [bedpe] ]
    versions    = ch_version                                // channel: [ path(version) ]
}
