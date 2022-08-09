/*
 * pair the proper mapped pairs
 */

include { COOLER_BALANCE } from '../../modules/local/cooler/balance/main'
include { COOLER_CLOAD   } from '../../modules/local/cooler/cload/main'
include { COOLER_MERGE   } from '../../modules/local/cooler/merge/main'
include { COOLER_ZOOMIFY } from '../../modules/local/cooler/zoomify/main'
include { COOLER_DUMP
    as COOLER_DUMP_PER_GROUP    } from '../../modules/nf-core/modules/cooler/dump/main'
include { COOLER_DUMP
    as COOLER_DUMP_PER_SAMPLE   } from '../../modules/nf-core/modules/cooler/dump/main'
include { DUMPREADS
    as DUMPREADS_PER_GROUP } from '../../modules/local/cooler/dumpreads'
include { DUMPREADS
    as DUMPREADS_PER_SAMPLE} from '../../modules/local/cooler/dumpreads'
include { JUICER         } from '../../modules/local/cooler/juicer'

workflow COOLER {
    take:
    valid_pairs               // channel: [ val(meta), val(bin), [pairs], [pairs.px] ]
    chromsizes                // channel: [ path(chromsizes) ]
    juicer_jvm_params         // values
    juicer_tools_jar          // channel: [ path(juicer_tool jar) ]
    long_bedpe_postfix

    main:
    // HiC-like contact matrix
    ch_version = COOLER_CLOAD(valid_pairs.map{[it[0], it[2], it[3]]}, valid_pairs.map{it[1]}, chromsizes).versions
    // Merge contacts
    COOLER_CLOAD.out.cool
                .map{
                    meta, bin, cool ->
                    [meta.group, bin, cool]
                }
                .groupTuple(by:[0, 1])
                .map{group, bin, cool -> [[id:group, bin:bin], cool]}
                .set{ch_cooler}
    COOLER_MERGE(ch_cooler)
    // create a balanced matrix for compartment and tad calls
    COOLER_BALANCE(COOLER_MERGE.out.cool)
    // create mcooler file for visualization
    COOLER_ZOOMIFY(COOLER_BALANCE.out.cool)
    // dump interaction bedpe for each group
    COOLER_DUMP_PER_GROUP(COOLER_MERGE.out.cool, [])
    // dump long interaction bedpe for each group for MAPS to call peaks
    DUMPREADS_PER_GROUP(COOLER_DUMP_PER_GROUP.out.bedpe, long_bedpe_postfix)
    if(juicer_tools_jar){
        JUICER(DUMPREADS_PER_GROUP.out.gi, juicer_tools_jar, chromsizes, juicer_jvm_params)
        ch_version = ch_version.mix(JUICER.out.versions)
    }

    // dump long interaction bedpe for each sample
    COOLER_DUMP_PER_SAMPLE(COOLER_CLOAD.out.cool.map{ meta, bin, cool -> [[id:meta.id, group:meta.group, bin:bin], cool]}, [])
    DUMPREADS_PER_SAMPLE(COOLER_DUMP_PER_SAMPLE.out.bedpe, long_bedpe_postfix)
    ch_version = ch_version.mix(DUMPREADS_PER_GROUP.out.versions)

    emit:
    cool        = COOLER_BALANCE.out.cool                   // channel: [ val(meta), [cool] ]
    mcool       = COOLER_ZOOMIFY.out.mcool                  // channel: [ val(meta), [mcool] ]
    groupbedpe  = COOLER_DUMP_PER_GROUP.out.bedpe           // channel: [ val(meta), [bedpe] ]
    bedpe       = DUMPREADS_PER_GROUP.out.bedpe             // channel: [ val(meta), [bedpe] ]
    samplebedpe = DUMPREADS_PER_SAMPLE.out.bedpe            // channel: [ val(meta), [bedpe] ]
    versions    = ch_version                                // channel: [ path(version) ]
}
