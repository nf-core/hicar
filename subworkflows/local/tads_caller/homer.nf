/*
 * Call TADs by Homer
 */

include { HOMER_FINDTADSANDLOOPS
    as HOMER_FINDTADSANDLOOPS_TADS        } from '../../../modules/local/homer/find_tad_loops'

workflow HOMER_TADS {
    take:
    tagdir       // channel: [ val(meta), [cool] ]
    bin_size     // channel: [ val(resolution) ]
    genome       // value: UCSC genome name

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()

    // call tads
    ch_tads = HOMER_FINDTADSANDLOOPS_TADS(
        tagdir.map{
            meta, tag -> [ meta, bin_size, tag ]
        },
        genome).tads
    ch_versions = HOMER_FINDTADSANDLOOPS_TADS.out.versions


    emit:
    tads         = ch_tads                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
