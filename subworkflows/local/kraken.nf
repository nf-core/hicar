//
// Detect the contamination
//

include { KRAKEN2_KRAKEN2             } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_COMBINEKREPORTS } from '../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_KREPORT2KRONA   } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT          } from '../../modules/nf-core/krona/ktimporttext/main'

workflow KRAKEN2 {
    take:
    reads4mapping           // reads [meta, [fastqs]]
    kraken2_db              // kraken2 database folder name [ [db] ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    KRAKEN2_KRAKEN2(reads4mapping, kraken2_db, false, false)

    ch_multiqc_files = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report.map{[it[1]]}.collect() )
    ch_versions = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.ifEmpty(null) )

    KRAKENTOOLS_COMBINEKREPORTS ( KRAKEN2_KRAKEN2.out.report.map{[it[1]]}.collect().map{[[id:"kraken2_combined_reports"], it]} )
    ch_versions = ch_versions.mix( KRAKENTOOLS_COMBINEKREPORTS.out.versions.ifEmpty(null) )

    KRAKENTOOLS_KREPORT2KRONA ( KRAKEN2_KRAKEN2.out.report )
    ch_versions = ch_versions.mix( KRAKENTOOLS_KREPORT2KRONA.out.versions.ifEmpty(null) )

    KRONA_KTIMPORTTEXT( KRAKENTOOLS_KREPORT2KRONA.out.txt.map{[it[1]]}.collect().map{[[id:"kraken2_combined_reports"], it]} )
    ch_versions = ch_versions.mix( KRONA_KTIMPORTTEXT.out.versions.ifEmpty(null) )

    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files     // channel: [ meta, report ]
}
