//
// Transcription factor enrichment analysis
//

include { HOMER_TFEA                 } from './tfea/homer'
include { BIOC_ATACSEQTFEA           } from '../../modules/local/bioc/atacseqtfea'

workflow TFEA {
    take:
    bed                     // peaks regions [meta, R1/2, [peak]]
    additional_param        // singals for each tools

    main:
    ch_versions             = Channel.empty()

    switch(params.tfea_tool){
        case "homer":
            HOMER_TFEA(bed, additional_param) // additional_param: genome
            ch_versions = ch_versions.mix(HOMER_TFEA.out.versions)
            break
        case "atacseqtfea":
            BIOC_ATACSEQTFEA(bed, additional_param)
            ch_versions = ch_versions.mix(BIOC_ATACSEQTFEA.out.versions)
            break
        default:
            HOMER_TFEA(bed, additional_param) // additional_param: genome
            ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME.out.versions)
            break
    }

    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
}
