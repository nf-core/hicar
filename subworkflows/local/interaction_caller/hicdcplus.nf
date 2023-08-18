/*
 * Call interaction peaks by HiC-CD+
 */

include { HICDCPLUS_FEATURES             } from '../../../modules/local/hicdcplus/features'
include { HICDCPLUS_CALL_LOOPS            } from '../../../modules/local/hicdcplus/callloops'

workflow HICDCPLUS {
    take:
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome reads for group
    reads                           // [ meta, [bedgraph] ] cooler dump reads for each group
    additional_param                // [
                                    //    0 values: site, eg. 'GATC 1'
                                    //    1 [ genome fa ],
                                    //    2 [ chromsizes ],
                                    //    3 [ mappability bigwig ],
                                    // ]

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()

    // generate features
    ch_versions = HICDCPLUS_FEATURES(
                    bedpe.map{[ it[0].bin ]} // bin_size
                        .unique() // unique the bin_size
                        .combine(additional_param)
                    ).versions
    // call loops
    ch_reads = reads.combine(HICDCPLUS_FEATURES.out.features)
                                    .map{[[id:it[0].id, bin:it[2]],
                                        it[1], it[3]]}
        .combine(bedpe, by: 0)
    ch_loop = HICDCPLUS_CALL_LOOPS(ch_reads.combine(additional_param.map{[it[2]]})).interactions
    ch_versions = ch_versions.mix(HICDCPLUS_CALL_LOOPS.out.versions.ifEmpty([]))

    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
