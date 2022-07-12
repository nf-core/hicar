/*
 * call TADs by HiCExploer
 */

include { HICEXPLORER_HICFINDTADS } from '../../modules/local/hicexplorer/hicfindtads'
include { HICEXPLORER_HICPLOTTADS } from '../../modules/local/hicexplorer/hicplottads'

workflow HICEXPLORER_CALLTADS {
    take:
    cool       // channel: [ val(meta), [cool] ]
    resolution // channel: [ val(resolution) ]
    chromsizes // channel: [ path(size) ]

    main:
    HICEXPLORER_HICFINDTADS(
        cool,
        resolution
    )
    ch_version = HICEXPLORER_HICFINDTADS.out.versions.ifEmpty(null)

    HICEXPLORER_HICPLOTTADS(
        cool.join(HICEXPLORER_HICFINDTADS.out.domains),
        chromsizes
    )
    ch_version = ch_version.mix(HICEXPLORER_HICPLOTTADS.out.versions.ifEmpty(null))

    emit:
    results   = HICEXPLORER_HICFINDTADS.out.results  // channel: [ val(meta), path(results)]
    domains   = HICEXPLORER_HICFINDTADS.out.domains  // channel: [ val(meta), path(domains.bed)]
    pngs      = HICEXPLORER_HICPLOTTADS.out.pngs     // channel: [ val(meta), path(pngs)]
    ini       = HICEXPLORER_HICPLOTTADS.out.ini      // channel: [ val(meta), path(ini)]
    versions  = ch_version                           // channel: [ path(version) ]
}
