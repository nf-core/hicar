/*
 * call TADs by HiCExplorer
 */

include { HICEXPLORER_HICFINDTADS } from '../../../modules/local/hicexplorer/hicfindtads'
include { HICEXPLORER_HICPLOTTADS } from '../../../modules/local/hicexplorer/hicplottads'

workflow HICEXPLORER_TADS {
    take:
    cool       // channel: [ val(meta), [cool] ]
    resolution // channel: [ val(resolution) ]
    chromsizes // channel: [ path(size) ]

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions      = Channel.empty()

    HICEXPLORER_HICFINDTADS(
        cool,
        resolution
    )
    ch_version = HICEXPLORER_HICFINDTADS.out.versions

    HICEXPLORER_HICPLOTTADS(
        cool.join(HICEXPLORER_HICFINDTADS.out.tads),
        chromsizes
    )
    ch_version = ch_version.mix(HICEXPLORER_HICPLOTTADS.out.versions)

    emit:
    tads      = HICEXPLORER_HICFINDTADS.out.tads     // channel: [ val(meta), val(bin), path(domains.bed)]
    mqc       = ch_multiqc_files                     // channel: [ path(mqc) ]
    versions  = ch_version                           // channel: [ path(version) ]
}
