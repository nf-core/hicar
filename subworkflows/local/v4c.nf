//
// virtual_4c
//

include { COOLTOOLS_V4C                    } from './v4c/cooltools'
include { HICEXPLORER_CHICPLOTVIEWPOINT    } from '../../modules/local/hicexplorer/chicplotviewpoint'
include { BIOC_TRACKVIEWER                 } from '../../modules/local/bioc/trackviewer'

workflow V4C {
    take:
    matrix                       // tuple val(bin_size), path(events), path(mcools), path(anchors.bed)
    raw_pairs                    // .unselected.pairs.gz of samfrag
    gtf
    chrom_sizes
    restrict

    main:
    ch_versions             = Channel.empty()

    // use anchor bed as reference points file
    switch(params.v4c_tool){
        case "cooltools":
            COOLTOOLS_V4C(
                matrix
            )
            ch_versions = ch_versions.mix(COOLTOOLS_V4C.out.versions)
            break
        case "hicexplorer":
            HICEXPLORER_CHICPLOTVIEWPOINT(
                matrix
            )
            ch_versions = ch_versions.mix(HICEXPLORER_CHICPLOTVIEWPOINT.out.versions)
            break
        case "trackviewer":
            BIOC_TRACKVIEWER(
                matrix,
                raw_pairs,
                gtf,
                chrom_sizes,
                restrict
            )
            ch_versions = ch_versions.mix(BIOC_TRACKVIEWER.out.versions)
            break
        default:
            BIOC_TRACKVIEWER(
                matrix,
                raw_pairs,
                gtf,
                chrom_sizes,
                restrict
            )
            ch_versions = ch_versions.mix(BIOC_TRACKVIEWER.out.versions)
            break
    }


    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
}
