//
// Call A/B compartment
//

include { COOLTOOLS_COMPARTMENTS } from '../../modules/local/cooltools/eigs-cis'
include { HOMCER_COMPARTMENTS    } from './compartments_caller/homer'

workflow COMPARTMENTS {
    take:
    matrix                                       // tuple val(meta), path(cool)
    tagdir                                       // tuple val(meta), path(tagdir)
    resolution
    fasta
    chromsizes
    genome

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_compartments         = Channel.empty() // a bigwig files channel

    switch(params.compartments_tool){
        case "cooltools":
            COOLTOOLS_COMPARTMENTS(
                matrix,
                resolution,
                fasta,
                chromsizes
            )
            ch_compartments = COOLTOOLS_COMPARTMENTS.out.compartments
            ch_versions = COOLTOOLS_COMPARTMENTS.out.versions
            ch_circos_files = COOLTOOLS_COMPARTMENTS.out.compartments
            break
        case "homer":
            HOMCER_COMPARTMENTS(
                tagdir,
                genome
            )
        default:
            COOLTOOLS_COMPARTMENTS(
                matrix,
                resolution,
                fasta,
                chromsizes
            )
            ch_compartments = COOLTOOLS_COMPARTMENTS.out.compartments
            ch_versions = COOLTOOLS_COMPARTMENTS.out.versions
            ch_circos_files = COOLTOOLS_COMPARTMENTS.out.compartments
            break
    }


    emit:
    compartments    = ch_compartments      // channel: [ meta, [compartments] ]
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
