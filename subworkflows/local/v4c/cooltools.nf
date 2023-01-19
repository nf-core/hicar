/*
 * virtual_4c by cooltools
 */

include { BIOC_SUBSETLOOPS    } from '../../../modules/local/bioc/subsetloops'
include { VIRTUAL4C_BY_COOLTOOLS } from '../../../modules/local/cooltools/virtual4c'

workflow COOLTOOLS_V4C {
    take:
    matrix        // [ bin_size, [cool], [de_loops], [viewpoint_pe]]

    main:
    ch_versions = BIOC_SUBSETLOOPS(
        matrix.map{it[3]},
        matrix.map{[it[0], it[2]]}).versions
    VIRTUAL4C_BY_COOLTOOLS(
        matrix.map{[it[0], it[1]]}
            .combine(BIOC_SUBSETLOOPS.out.bed, by: 0),
        params.v4c_max_events
    )
    ch_versions = ch_versions.mix(VIRTUAL4C_BY_COOLTOOLS.out.versions)

    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
    v4c             = VIRTUAL4C_BY_COOLTOOLS.out.v4c
}
