/*
 * call TADs by HiCExplorer
 */

include { HICEXPLORER_HICPCA          } from '../../../modules/local/hicexplorer/hicpca'
include { HICEXPLORER_HICTRANSFORM    } from '../../../modules/local/hicexplorer/hictransform'
include { HICEXPLORER_HICPLOTMATRIX
    as HICEXPLORER_HICPLOTMATRIX_PCA1 } from '../../../modules/local/hicexplorer/hicplotmatrix'
include { HICEXPLORER_HICPLOTMATRIX
    as HICEXPLORER_HICPLOTMATRIX_PCA2 } from '../../../modules/local/hicexplorer/hicplotmatrix'

workflow HICEXPLORER_COMPARTMENTS {
    take:
    cool       // channel: [ val(meta), [cool] ]
    resolution // channel: [ val(resolution) ]
    chromsizes // channel: [ path(size) ]

    main:
    ch_versions      = Channel.empty()

    // step1 create bigwig files, the input file is balanced by cooler::balance (aka, ICE)
    HICEXPLORER_HICPCA(
        cool
    )
    ch_version = HICEXPLORER_HICPCA.out.versions

    // step2 create intermeidated files for plot
    HICEXPLORER_HICTRANSFORM(
        cool
    )
    ch_version = ch_version.mix(HICEXPLORER_HICTRANSFORM.out.versions)

    // step3 plot the compartments for pca1 and pca2
    HICEXPLORER_HICPLOTMATRIX_PCA1(
        HICEXPLORER_HICTRANSFORM.out.transformed
            .combine(HICEXPLORER_HICPCA.out.pca1, by: 0)
    )
    HICEXPLORER_HICPLOTMATRIX_PCA2(
        HICEXPLORER_HICTRANSFORM.out.transformed
            .combine(HICEXPLORER_HICPCA.out.pca2, by: 0)
    )
    ch_version = ch_version.mix(HICEXPLORER_HICPLOTMATRIX_PCA1.out.versions)

    emit:
    compartments  = HICEXPLORER_HICPCA.out.pca          // channel: [ val(meta), path(bigwig)]
    versions      = ch_version                           // channel: [ path(version) ]
}
