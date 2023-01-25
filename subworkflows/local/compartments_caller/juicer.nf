/*
 * call compartments by Juicer
 */

include { JUICER_EIGENVECTOR      } from '../../../modules/local/juicer/eigenvector'

workflow JUICER_COMPARTMENTS {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    resolution             // value
    additional_param       // value: [ [juicer_box_jar], [chrom_size] ]

    main:
    JUICER_EIGENVECTOR(
        matrix.combine(additional_param),
        resolution
    )
    ch_version = JUICER_EIGENVECTOR.out.versions


    emit:
    compartments   = JUICER_EIGENVECTOR.out.compartments         // channel: [ val(meta), path(bigwig)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
