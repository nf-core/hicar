/*
 * Call interaction peaks by HiCExplorer, current not proper
 */

include { HICEXPLORER_HICDETECTLOOPS      } from '../../../modules/local/hicexplorer/hicdetectloops'

workflow HICEXPLORER_INTERACTIONS {
    take:
    matrix                          // [ meta, [matrix] ] corrected matrix file

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()
    ch_versions = Channel.empty()

    // call interactions
    ch_loop = HICEXPLORER_HICDETECTLOOPS(matrix).bedpe
    ch_versions = HICEXPLORER_HICDETECTLOOPS.out.versions


    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
