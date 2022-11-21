//
// Call interaction loops
// input:
// [[2D signals], [1D peaks], [additional_files]]
// Ask all the output of each module must contain channel:
// versions: [path]; interactions: [ meta, bin_size, [bedpe] ]
// optional output: mqc: [path];
//

include { HICDCPLUS          } from './interaction_caller/hicdcplus'
include { MAPS               } from './interaction_caller/maps'
include { MERGE_INTERACTIONS } from '../../modules/local/bioc/merge_interactions'

workflow INTERACTIONS {
    take:
    matrix
    peaks
    additional_param

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_loops                = Channel.empty() // a bed files channel

    switch(params.interactions_tool){
        case "maps":
            MAPS(
                matrix,
                peaks,
                additional_param
            )
            ch_loops = MAPS.out.interactions
            ch_versions = MAPS.out.versions
            ch_annotation_files = MAPS.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = MAPS.out.interactions.map{
                meta, bin_size, bedpe ->
                    meta.bin = bin_size
                    [meta, bedpe]
            }
            break
        case "hicdcplus":
            HICDCPLUS(
                matrix,
                peaks,
                additional_param
            )
            ch_loops = HICDCPLUS.out.interactions
            ch_versions = HICDCPLUS.out.versions
            ch_annotation_files = HICDCPLUS.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = HICDCPLUS.out.interactions.map{
                meta, bin_size, bedpe ->
                    meta.bin = bin_size
                    [meta, bedpe]
            }
            break
    }
    // merge loops
    MERGE_INTERACTIONS(
        ch_loops.map{
            meta, bin_size, interactions -> [bin_size, interactions]
        }.groupTuple()
    )
    ch_versions = ch_versions.mix(MERGE_INTERACTIONS.out.versions.ifEmpty(null))

    emit:
    loops           = ch_loops                              // channel: [ meta, bin_size, path(bedpe) ]
    mergedloops     = MERGE_INTERACTIONS.out.interactions   // channel: [ bin_siz, path(bedpe) ]
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions                           // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
