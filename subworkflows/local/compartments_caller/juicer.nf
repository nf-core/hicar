/*
 * call compartments by Juicer
 */

include { JUICER_EIGENVECTOR      } from '../../../modules/local/juicer/eigenvector'
include { UCSC_WIGTOBIGWIG
    as COMPARTMENTS_WIGTOBIGWIG   } from '../../../modules/nf-core/ucsc/wigtobigwig/main'


workflow JUICER_COMPARTMENTS {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    resolution             // value
    additional_param       // value: [ [juicer_tools_jar], [chrom_size] ]

    main:
    JUICER_EIGENVECTOR(
        matrix.combine(additional_param),
        resolution
    )
    ch_version = JUICER_EIGENVECTOR.out.versions

    COMPARTMENTS_WIGTOBIGWIG(
        JUICER_EIGENVECTOR.out.wig,
        additional_param.map{it[1]}.collect()
    )
    ch_version = ch_version.mix(COMPARTMENTS_WIGTOBIGWIG.out.versions)

    emit:
    eigenvectors   = JUICER_EIGENVECTOR.out.compartments         // channel: [ val(meta), path(bigwig)]
    compartments   = COMPARTMENTS_WIGTOBIGWIG.out.bw             // channel: [ val(meta), path(bigwig)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
