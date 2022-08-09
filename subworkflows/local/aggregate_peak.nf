//
// aggregate peak analysis (APA)
//

include { HICEXPLORER_HICAGGREGATECONTACTS } from '../../modules/local/hicexplorer/hicaggregatecontacts'

workflow APA {
    take:
    matrix                                       // tuple val(meta), path(cool)
    peaks

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_apa                  = Channel.empty() // a png files channel

    switch(params.apa_tool){
        case "hicaggregatecontacts":
            HICEXPLORER_HICAGGREGATECONTACTS(
                matrix,
                peaks
            )
            ch_apa = HICEXPLORER_HICAGGREGATECONTACTS.out.png
            ch_versions = HICEXPLORER_HICAGGREGATECONTACTS.out.versions
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
