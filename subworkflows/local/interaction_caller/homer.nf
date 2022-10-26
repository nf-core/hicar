/*
 * Call interaction peaks by Homer, current not proper
 */

include { HOMER_ANALYZEHIC             } from '../../../modules/local/homer/analyze_hic'

workflow HOMER_INTERACTIONS {
    take:
    tagdir                          // [ meta, [tagdir] ] cooler dump ATAC reads for each group
    bin_size                        // values: bin size, 5000, 10000

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()
    ch_versions = Channel.empty()

    // call interactions
    tag = tagdir.combine(bin_size).map{
        meta, tag, bin -> [ meta, bin, tag]
    }
    ch_loop = HOMER_ANALYZEHIC(tag).bedpe
    ch_versions = HOMER_ANALYZEHIC.out.versions


    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
