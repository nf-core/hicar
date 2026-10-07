/*
 * call APA by Juicer
 */

include { JUICER_APA       } from '../../../modules/local/juicer/apa'
include { BIOC_SUBSETLOOPS } from '../../../modules/local/bioc/subsetloops'

workflow JUICER_APACALLER {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    peaks                  // path(1D peaks)
    additional_param       // values [bin, merged_loops, juicer_tools_jar]
    format                 // output plot format

    main:
    ch_mergedloops = additional_param.map{[it[0], it[1]]}
    BIOC_SUBSETLOOPS(peaks, ch_mergedloops)

    ch_hic_loops = matrix.map{[it[0].bin, it[0], it[1]]}
                        .combine(BIOC_SUBSETLOOPS.out.loops, by: 0)
                        .map{[it[1], it[2], it[3]]} // merge by bin

    JUICER_APA(
        ch_hic_loops,
        additional_param.map{[it[2]]}
    )
    ch_version = JUICER_APA.out.versions


    emit:
    plot      = JUICER_APA.out.plot                              // channel: [ val(meta), path(pngs)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
