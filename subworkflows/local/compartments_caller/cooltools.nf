/*
 * call compartments by cooltools
 */

include { COOLTOOLS_EIGSCIS   } from '../../../modules/local/cooltools/eigscis'

workflow COOLTOOLS_COMPARTMENTS {
    take:
    maxtix                 // channel: [ val(meta), [cool] ]
    resolution             // value
    additional_param       // values [fasta, chromosizes]

    main:
    COOLTOOLS_EIGSCIS(
        maxtix.combine(additional_param),
        resolution
    )
    ch_version = COOLTOOLS_EIGSCIS.out.versions


    emit:
    compartments   = COOLTOOLS_EIGSCIS.out.compartments          // channel: [ val(meta), path(bigwig)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
