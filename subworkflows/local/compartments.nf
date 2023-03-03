//
// Call A/B compartment
//

include { COOLTOOLS_COMPARTMENTS   } from './compartments_caller/cooltools'
include { HOMER_COMPARTMENTS       } from './compartments_caller/homer'
include { JUICER_COMPARTMENTS      } from './compartments_caller/juicer'
include { HICEXPLORER_COMPARTMENTS } from './compartments_caller/hicexplorer'

workflow COMPARTMENTS {
    take:
    matrix                                       // tuple val(meta), path(cool/tagdir)
    resolution
    additional_param
    cool_bin

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty() // TODO
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_compartments         = Channel.empty() // a compartments files channel

    switch(params.compartments_tool){
        case "cooltools":
            COOLTOOLS_COMPARTMENTS(
                matrix,
                resolution,
                additional_param //fasta, chrom_size
            )
            ch_compartments = COOLTOOLS_COMPARTMENTS.out.compartments
            ch_versions = COOLTOOLS_COMPARTMENTS.out.versions
            ch_circos_files = COOLTOOLS_COMPARTMENTS.out.compartments
            break
        case "homer":
            HOMER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param  //genome, chrom_size
            )
            ch_compartments = HOMER_COMPARTMENTS.out.compartments
            ch_versions = HOMER_COMPARTMENTS.out.versions
            ch_circos_files = HOMER_COMPARTMENTS.out.compartments.combine(cool_bin).map{[[id:it[0].id, bin:it[2]], it[1]]}
            break
        case "juicebox":
            JUICER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param
            )
            ch_compartments = JUICER_COMPARTMENTS.out.compartments
            ch_versions = JUICER_COMPARTMENTS.out.versions
            ch_circos_files = JUICER_COMPARTMENTS.out.compartments
            break
        case "hicexplorer":
            HICEXPLORER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param // chromsizes
            )
            ch_tads = HICEXPLORER_COMPARTMENTS.out.compartments
            ch_versions = HICEXPLORER_COMPARTMENTS.out.versions
            ch_circos_files = HICEXPLORER_COMPARTMENTS.out.compartments
            break
        default:
            COOLTOOLS_COMPARTMENTS(
                matrix,
                resolution,
                additional_param //fasta, chrom_size
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
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
