//
// aggregate peak analysis (APA)
// The input is matrix and peaks (1D/2D)
// Ask all the output of each module must contain channel:
// versions: [path]; png: [ meta, [png] ]
// optional output: mqc: [path];
//

include { HICEXPLORER_HICAGGREGATECONTACTS } from '../../modules/local/hicexplorer/hicaggregatecontacts'
include { JUICER_APA                       } from '../../modules/local/juicer/apa'

workflow APA {
    take:
    matrix            // tuple val(meta), path(cool/hic), signal file
    peaks             // path(1d/2d peaks) for APA analysis

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_apa                  = Channel.empty() // a png files channel

    switch(params.apa_tool){
        case "hicexploer":
            HICEXPLORER_HICAGGREGATECONTACTS(
                matrix,
                peaks
            )
            ch_apa = HICEXPLORER_HICAGGREGATECONTACTS.out.png
            ch_versions = HICEXPLORER_HICAGGREGATECONTACTS.out.versions
            break
        case "juicebox":
            JUICER_APA(
                matrix,
                peaks
            )
            ch_apa = JUICER_APA.out.png
            ch_versions = JUICER_APA.out.versions
            break
        default:
            HICEXPLORER_HICAGGREGATECONTACTS(
                matrix,
                peaks
            )
            ch_apa = HICEXPLORER_HICAGGREGATECONTACTS.out.png
            ch_versions = HICEXPLORER_HICAGGREGATECONTACTS.out.versions
            break
    }


    emit:
    apa             = ch_apa               // channel: [ meta, [png] ]
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
