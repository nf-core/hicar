//
// Call TADs
//

include { COOLTOOLS_INSULATION } from '../../modules/local/cooltools/insulation'
include { HICEXPLORER_CALLTADS } from './tads_caller/hicexplorer'

workflow TADS {
    take:
    matrix                                       // tuple val(meta), path(cool)
    resolution
    fasta
    chromsizes

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_tads                 = Channel.empty() // a bed files channel

    switch(params.tad_tool){
        case "insulation":
            COOLTOOLS_INSULATION(
                matrix,
                resolution
            )
            ch_tads = COOLTOOLS_INSULATION.out.tads
            ch_versions = COOLTOOLS_INSULATION.out.versions
            ch_circos_files = COOLTOOLS_INSULATION.out.tads
            break
        case "hicfindtads":
            HICEXPLORER_CALLTADS(
                matrix,
                resolution,
                chromsizes
            )
            ch_tads = HICEXPLORER_CALLTADS.out.tads
            ch_versions = HICEXPLORER_CALLTADS.out.versions
            ch_circos_files = HICEXPLORER_CALLTADS.out.tads
            break
        default:
            HICEXPLORER_CALLTADS(
                matrix,
                resolution,
                chromsizes
            )
            ch_tads = HICEXPLORER_CALLTADS.out.tads
            ch_versions = HICEXPLORER_CALLTADS.out.versions
            ch_circos_files = COOLTOOLS_INSULATION.out.tads
            break
    }


    emit:
    tads            = ch_tads              // channel: [ meta, [TADs] ]
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
