/*
 * virtual_4c by cooltools
 */

include { BIOC_SUBSETLOOPS    } from '../../../modules/local/bioc/subsetloops'
include { COOLTOOLS_VIRTUAL4C } from '../../../modules/local/cooltools/virtual4c'

workflow COOLTOOLS_V4C {
    take:
    matrix        // [ bin_size, [cool], [de_loops], [viewpoint_pe]]

    main:
    ch_versions = BIOC_SUBSETLOOPS(
        matrix.map{it[3]},
        matrix.map{[it[0], it[2]]}).versions
    COOLTOOLS_VIRTUAL4C(
        matrix.map{[it[0], it[1]]}
            .combine(BIOC_SUBSETLOOPS.out.bed, by: 0),
        params.v4c_max_events
    )
    ch_versions = ch_versions.mix(COOLTOOLS_VIRTUAL4C.out.versions)

    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
    v4c             = COOLTOOLS_VIRTUAL4C.out.v4c
}
