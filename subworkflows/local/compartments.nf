//
// Call A/B compartment
//

include { COOLTOOLS_COMPARTMENTS   } from './compartments_caller/cooltools'
include { HOMER_COMPARTMENTS       } from './compartments_caller/homer'
include { JUICER_COMPARTMENTS      } from './compartments_caller/juicer'
include { HICEXPLORER_COMPARTMENTS } from './compartments_caller/hicexplorer'
include { ADJUST_COMPARTMENTS      } from '../../modules/local/bioc/adjust_compartments'
include { DIFFERENTIAL_COMPARTMENTS} from '../../modules/local/bioc/differential_compartments'

workflow COMPARTMENTS {
    take:
    matrix                                       // tuple val(meta), path(cool/tagdir)
    resolution
    additional_param
    cool_bin
    //atac_peaks                                   // path(bed) ATAC accessibility could be used to adjust A/B compartments

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
            break
        case "homer":
            HOMER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param  //genome, chrom_size
            )
            ch_compartments = HOMER_COMPARTMENTS.out.compartments.combine(cool_bin).map{[[id:it[0].id, bin:it[2]], it[1]]}
            ch_versions = HOMER_COMPARTMENTS.out.versions
            break
        case "juicebox":
            JUICER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param
            )
            ch_compartments = JUICER_COMPARTMENTS.out.compartments
            ch_versions = JUICER_COMPARTMENTS.out.versions
            break
        case "hicexplorer":
            HICEXPLORER_COMPARTMENTS(
                matrix,
                resolution,
                additional_param // chromsizes
            )
            ch_compartments = HICEXPLORER_COMPARTMENTS.out.compartments
            ch_versions = HICEXPLORER_COMPARTMENTS.out.versions
            break
        default:
            COOLTOOLS_COMPARTMENTS(
                matrix,
                resolution,
                additional_param //fasta, chrom_size
            )
            ch_compartments = COOLTOOLS_COMPARTMENTS.out.compartments
            ch_versions = COOLTOOLS_COMPARTMENTS.out.versions
            break
    }
    // adjust A/B compartments
    ADJUST_COMPARTMENTS(ch_compartments.map{meta,bw -> [meta.bin, meta, bw]}.groupTuple())
    ch_compartments = ADJUST_COMPARTMENTS.out.compartments.transpose().map{bin, meta, bw -> [meta, bw]}
    ch_versions = ch_versions.mix(ADJUST_COMPARTMENTS.out.versions)
    ch_circos_files = ADJUST_COMPARTMENTS.out.compartments.transpose().map{bin, meta, bw -> [meta, bw]}
    ch_track_files = ADJUST_COMPARTMENTS.out.compartments.transpose().map{bin, meta, bw -> [meta, bw]}

    // differential analysis
    DIFFERENTIAL_COMPARTMENTS(ADJUST_COMPARTMENTS.out.compartments)
    ch_versions = DIFFERENTIAL_COMPARTMENTS.out.versions
    ch_track_files = ch_track_files.mix(DIFFERENTIAL_COMPARTMENTS.out.diff.map{
        bin, bw -> [[id: "differential_compartments_$bin"], bw]})

    emit:
    compartments    = ch_compartments      // channel: [ meta, [compartments] ]
    circos          = ch_circos_files      // channel: [ meta, [compartments] ]
    igv             = ch_track_files       // channel: [ meta, [compartments] ]
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
