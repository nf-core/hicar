/*
 * TFEA by homer
 */

include { HOMER_FINDMOTIFSGENOME     } from '../../../modules/local/homer/find_motifs_genome'
include { SEQLEVELS_STYLE            } from '../../../modules/local/bioc/seqlevelsstyle'
include { ENSEMBL_UCSC_CONVERT       } from '../../../modules/local/bioc/ensembl_ucsc_convert'

workflow HOMER_TFEA {
    take:
    bed                      // peaks regions [meta, R1/2, [peak]]
    additional_param         // signals for each tools

    main:
    ch_versions = Channel.empty()

    if(!workflow.containerEngine){
        // check seqlevels
        seqlevelsstyle = SEQLEVELS_STYLE(bed.map{it[2]}.collect().map{it[0]}).seqlevels_style
        ch_versions = ch_versions.mix(SEQLEVELS_STYLE.out.versions)
        if("$seqlevelsstyle" != "UCSC"){
            ch_bed = bed.map{[it[0].id, it[0], it[1], it[2]]}
            ch_new_bed = ch_bed.combine(ENSEMBL_UCSC_CONVERT(ch_bed.map{[it[0], it[3]]}).tab, by:0)
                            .map{[it[1], it[2], it[4]]}
            ch_versions = ch_versions.mix(ENSEMBL_UCSC_CONVERT.out.versions)
        }
    }else{
        ch_new_bed = bed
    }
    HOMER_FINDMOTIFSGENOME(ch_new_bed, additional_param)
    ch_versions = ch_versions.mix(HOMER_FINDMOTIFSGENOME.out.versions)


    emit:
    versions     = ch_versions                  // channel: [ path(version) ]
}
