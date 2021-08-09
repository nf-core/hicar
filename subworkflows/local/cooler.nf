/*
 * pair the proper mapped pairs
 */
params.options = [:]

include { COOLER_CLOAD   } from '../../modules/local/cooler/cload/main'   addParams(options: params.options.cooler_cload)
include { COOLER_MERGE   } from '../../modules/local/cooler/merge/main'   addParams(options: params.options.cooler_merge)
include { COOLER_ZOOMIFY } from '../../modules/local/cooler/zoomify/main' addParams(options: params.options.cooler_zoomify)
include { COOLER_DUMP    } from '../../modules/nf-core/modules/cooler/dump/main'    addParams(options: params.options.cooler_dump_per_group)
include { COOLER_DUMP
    as COOLER_DUMP_SAMPLE} from '../../modules/nf-core/modules/cooler/dump/main'    addParams(options: params.options.cooler_dump_per_sample)
include { DUMPINTRAREADS } from '../../modules/local/cooler/dumpintrareads'    addParams(options: params.options.dumpintrareads_per_group)
include { DUMPINTRAREADS
    as DUMPINTRAREADS_SAMPLE} from '../../modules/local/cooler/dumpintrareads'    addParams(options: params.options.dumpintrareads_per_sample)

workflow COOLER {
    take:
    valid_pairs  // channel: [ val(meta), val(bin), [pairs], [pairs.px] ]
    chromsizes   // channel: [ path(chromsizes) ]

    main:
    // HiC-like contact matrix
    ch_version = COOLER_CLOAD(valid_pairs.map{[it[0], it[2], it[3]]}, valid_pairs.map{it[1]}, chromsizes).version
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
    // create mcooler file for visualization
    COOLER_ZOOMIFY(COOLER_MERGE.out.cool)
    // dump long.intra.bedpe for each group for MAPS to call peaks
    COOLER_DUMP(COOLER_MERGE.out.cool).bedpe | DUMPINTRAREADS
    // dump long.intra.bedpe for each sample
    COOLER_DUMP_SAMPLE(COOLER_CLOAD.out.cool.map{ meta, bin, cool -> [[id:meta.id, group:meta.group, bin:bin], cool]})
    DUMPINTRAREADS_SAMPLE(COOLER_DUMP_SAMPLE.out.bedpe)
    ch_version = ch_version.mix(DUMPINTRAREADS.out.version)

    emit:
    mcool       = COOLER_ZOOMIFY.out.mcool        // channel: [ val(meta), [mcool] ]
    bedpe       = DUMPINTRAREADS.out.bedpe        // channel: [ val(meta), [bedpe] ]
    samplebedpe = DUMPINTRAREADS_SAMPLE.out.bedpe // channel: [ val(meta), [bedpe] ]
    version     = ch_version                      // channel: [ path(version) ]
}
