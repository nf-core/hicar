//
// Get coverage scale
//

include { GET_SCALE           } from '../../modules/local/atacreads/getscale'

workflow COVERAGE_SCALE {
    take:
    counts                     // counts [meta, counts]

    main:
    ch_versions             = Channel.empty()

    min_cnt = counts.map{[it[1].toInteger()]}.collect().map{it.min()}
    GET_SCALE(counts, min_cnt)

    emit:
    scale           = GET_SCALE.out.scale
    versions        = ch_versions          // channel: [ versions.yml ]
}
