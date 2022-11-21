/*
 * call APA by Juicer
 */

include { JUICER_APA       } from '../../../modules/local/juicer/apa'
include { BIOC_SUBSETLOOPS } from '../../../modules/local/bioc/subsetloops'

workflow JUICER_APACALLER {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    peaks                  // path(1D peaks)
    additional_param       // values [juicer_box_jar, chromsize, bin, merged_loops]

    main:
    //additional_param.view()
    ch_mergedloops = additional_param.map{it[3]}
    //ch_mergedloops.view()
    BIOC_SUBSETLOOPS(peaks, ch_mergedloops)

    ch_hic_loops = matrix.map{[it[0].bin, it[0], it[1]]}
                        .combine(BIOC_SUBSETLOOPS.out.loops, by: 0)
                        .map{[it[1], it[2], it[3]]} // merge by bin

    //ch_hic_loops.view{"loops: $it"}
    //peaks.view{"peak: $it"}
    //additional_param.map{[it[0], it[1]]}.view()
    JUICER_APA(
        ch_hic_loops,
        additional_param.map{[it[0]]},
        params.juicer_jvm_params
    )
    ch_version = JUICER_APA.out.versions.ifEmpty(null)


    emit:
    png       = JUICER_APA.out.png                               // channel: [ val(meta), path(pngs)]
    versions  = ch_version                                       // channel: [ path(version) ]
}