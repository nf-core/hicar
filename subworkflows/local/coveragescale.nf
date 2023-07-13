//
// Get coverage scale
//

workflow COVERAGE_SCALE {
    take:
    counts                     // counts [meta, counts]

    main:
    ch_versions             = Channel.empty()

    min_cnt = counts.map{[it[1].toInteger()]}.collect().map{it.min()}
    scale = counts.combine(min_cnt).map{
        meta, counts, min_counts ->
            factor = min_counts/counts.toInteger()
            [meta, factor]
    }
    scale.view()

    emit:
    scale                                  // channel: [ meta, scale_factor ]
    versions        = ch_versions          // channel: [ versions.yml ]
}
