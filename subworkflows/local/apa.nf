//
// aggregate peak analysis (APA)
// The input is matrix and peaks (1D/2D)
// Ask all the output of each module must contain channel:
// versions: [path]; png: [ meta, [png] ]
// optional output: mqc: [path];
//

include { HICEXPLORER_HICAGGREGATECONTACTS } from '../../modules/local/hicexplorer/hicaggregatecontacts'
include { JUICER_APACALLER                 } from './apa/juicer'
include { COOLTOOLS_APACALLER              } from './apa/cooltools'

workflow APA {
    take:
    matrix            // tuple val(meta), path(cool/hic), signal file
    peaks             // path(1d/2d peaks) for APA analysis
    additional_param

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty() // TODO
    ch_apa                  = Channel.empty() // a png files channel

    switch(params.apa_tool){
        case "hicexplorer":
            HICEXPLORER_HICAGGREGATECONTACTS(
                matrix,
                peaks,
                params.apa_format
            )
            ch_apa = HICEXPLORER_HICAGGREGATECONTACTS.out.plot
            ch_versions = HICEXPLORER_HICAGGREGATECONTACTS.out.versions
            break
        case "juicebox":
            JUICER_APACALLER(
                matrix,
                peaks,
                additional_param,
                params.apa_format
            )
            ch_apa = JUICER_APACALLER.out.plot
            ch_versions = JUICER_APACALLER.out.versions
            break
        case "cooltools":
            COOLTOOLS_APACALLER(
                matrix,
                peaks,
                params.apa_format
            )
            ch_apa = COOLTOOLS_APACALLER.out.plot
            ch_versions = COOLTOOLS_APACALLER.out.versions
            break
        default:
            HICEXPLORER_HICAGGREGATECONTACTS(
                matrix,
                peaks,
                params.apa_format
            )
            ch_apa = HICEXPLORER_HICAGGREGATECONTACTS.out.plot
            ch_versions = HICEXPLORER_HICAGGREGATECONTACTS.out.versions
            break
    }


    emit:
    apa             = ch_apa               // channel: [ meta, [plot] ]
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
