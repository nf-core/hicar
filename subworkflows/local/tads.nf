//
// Call TADs
// Ask all the output of each module must contain channel:
// versions: [path]; tads: [ meta, bin_size, [bedpe] ]
// optional output: mqc: [path];
//

include { COOLTOOLS_TADS       } from './tads_caller/cooltools'
include { HICEXPLORER_TADS     } from './tads_caller/hicexplorer'
include { HOMER_TADS           } from './tads_caller/homer'

workflow TADS {
    take:
    matrix                                       // tuple val(meta), path(cool/tagdir)
    resolution                                   // resolution for TADs calling
    additional_param                             // additional parameters

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty() // TODO
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_tads                 = Channel.empty() // a bed files channel

    switch(params.tad_tool){
        case "cooltools":
            COOLTOOLS_TADS(
                matrix,
                resolution
            )
            ch_tads = COOLTOOLS_TADS.out.tads
            ch_versions = COOLTOOLS_TADS.out.versions
            ch_circos_files = COOLTOOLS_TADS.out.tads.map{[it[0], it[2]]}
            break
        case "hicexploer":
            HICEXPLORER_TADS(
                matrix,
                resolution,
                additional_param  //chromsizes
            )
            ch_tads = HICEXPLORER_TADS.out.tads
            ch_versions = HICEXPLORER_TADS.out.versions
            ch_circos_files = HICEXPLORER_TADS.out.tads.map{[it[0], it[2]]}
            break
        case "homer":
            HOMER_TADS(
                matrix,
                resolution,
                additional_param //genome, ucsc genome name
            )
            ch_tads = HOMER_TADS.out.tads
            ch_versions = HOMER_TADS.out.versions
            ch_circos_files = HOMER_TADS.out.tads.map{[it[0], it[2]]}
            break
        default:
            HICEXPLORER_TADS(
                matrix,
                resolution,
                additional_param //chromsizes
            )
            ch_tads = HICEXPLORER_TADS.out.tads
            ch_versions = HICEXPLORER_TADS.out.versions
            ch_circos_files = HICEXPLORER_TADS.out.tads.map{[it[0], it[2]]}
            break
    }


    emit:
    tads            = ch_tads              // channel: [ meta, bin_size, [TADs] ]
    circos          = ch_circos_files      // channel: [ meta, [TADs]]
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
