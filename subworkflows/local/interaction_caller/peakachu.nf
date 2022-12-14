/*
 * Call interaction peaks by HiC-CD+
 */

include { PEAKACHU_MODEL             } from '../../../modules/local/peakachu/pretrained'
include { PEAKACHU_SCORE             } from '../../../modules/local/peakachu/score_genome'

workflow PEAKACHU {
    take:
    cool                           // [ val(meta), [cool] ] cool file

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()

    // generate predefined model file
    ch_versions = PEAKACHU_MODEL( cool ).versions
    // call loops
    ch_loop = PEAKACHU_SCORE(cool.combine(PEAKACHU_MODEL.out.model, by:0)).interactions
    ch_versions = ch_versions.mix(PEAKACHU_SCORE.out.versions.ifEmpty([]))

    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
