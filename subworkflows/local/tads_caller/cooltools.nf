/*
 * Call TADs by CoolTools
 */

 include { COOLTOOLS_INSULATION } from '../../../modules/local/cooltools/insulation'

workflow COOLTOOLS_TADS {
    take:
    matrix         // channel: [ val(meta), [cool] ]
    resolution     // channel: [ val(resolution) ]

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()

    // call tads
    ch_tads = COOLTOOLS_INSULATION(
        matrix,
        resolution
    ).tads
    ch_versions = COOLTOOLS_INSULATION.out.versions


    emit:
    tads         = ch_tads                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
