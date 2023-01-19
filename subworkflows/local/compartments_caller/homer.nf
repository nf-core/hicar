/*
 * call compartments by homer
 */

include { HOMER_RUNHICPCA } from '../../../modules/local/homer/run_hic_pca'
include { HOMER_FINDHICCOMPARTMENTS } from '../../../modules/local/homer/find_hic_compartments'

workflow HOMER_COMPARTMENTS {
    take:
    tagdir       // channel: [ val(meta), [tagdir] ]
    resolution   // value
    genome       // value

    main:
    HOMER_RUNHICPCA(
        tagdir,
        resolution,
        genome
    )
    ch_version = HOMER_RUNHICPCA.out.versions

    HOMER_FINDHICCOMPARTMENTS(
        HOMER_RUNHICPCA.out.txt,
        genome
    )
    ch_version = ch_version.mix(HOMER_FINDHICCOMPARTMENTS.out.versions)

    emit:
    compartments   = HOMER_FINDHICCOMPARTMENTS.out.compartments  // channel: [ val(meta), path(compartments)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
