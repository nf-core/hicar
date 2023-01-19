/*
 * call compartments by Juicer
 */

include { JUICER_EIGENVECTOR      } from '../../../modules/local/juicer/eigenvector'

workflow JUICER_COMPARTMENTS {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    resolution             // value
    additional_param       // values [jvm_params, juicer_box_jar]

    main:
    JUICER_EIGENVECTOR(
        matrix,
        resolution,
        additional_param
    )
    ch_version = JUICER_EIGENVECTOR.out.versions.ifEmpty(null)


    emit:
    compartments   = JUICER_EIGENVECTOR.out.compartments         // channel: [ val(meta), path(bigwig)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
