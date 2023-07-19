/*
 * Differentail analysis by HiCExplorer::Capture HiC pipeline
 */

include { HICEXPLORER_CHICSIGNIFICANTINTERACTIONS } from '../../../modules/local/hicexplorer/chicsignificantinteractions'
include { HICEXPLORER_CHICAGGREGATESTATISTIC } from '../../../modules/local/hicexplorer/chicaggregatestatistic'
include { HICEXPLORER_CHICDIFFERENTIALTEST } from '../../../modules/local/hicexplorer/chicdifferentialtest'

workflow HICEXPLORER_DIFFHIC {
    take:
    matrix        // [ bin_size, [background], [interactions] ]

    main:
    // significant target sites detection
    HICEXPLORER_CHICSIGNIFICANTINTERACTIONS(matrix)
    ch_versions = HICEXPLORER_CHICSIGNIFICANTINTERACTIONS.out.versions
    // aggregate data for differential test
    HICEXPLORER_CHICAGGREGATESTATISTIC(
        HICEXPLORER_CHICSIGNIFICANTINTERACTIONS.out.interaction_target
    )
    // differential test
    HICEXPLORER_CHICDIFFERENTIALTEST(
        HICEXPLORER_CHICAGGREGATESTATISTIC.out.aggregate
    )
    ch_diff = matrix.combine(HICEXPLORER_CHICDIFFERENTIALTEST.out.differential, by: 0)

    // differential analsis ask output anno [binsize, prefix, [files]], diff, stats, versions
    // the current hicExplorer does not export bedpe results, skip annotation.
    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
    anno            = Channel.empty()      // [bin_size, prefix, [files]]
    diff            = ch_diff              // [bin_size, [background], [interactions], [differential]] for v4c plot
}
