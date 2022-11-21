/*
 * Call interaction peaks by HiCExplorer, current not proper
 */

include { HICEXPLORER_CHICSIGNIFICANTINTERACTIONS      } from '../../../modules/local/hicexplorer/chicsignificantinteractions'


workflow HICEXPLORER_INTERACTIONS {
    take:
    matrix                          // [ bin_size, [cool] ] corrected matrix file

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()
    ch_versions = Channel.empty()

    ch_versions = HICEXPLORER_CHICSIGNIFICANTINTERACTIONS(matrix).versions

    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
