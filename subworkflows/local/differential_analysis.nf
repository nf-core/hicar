//
// differential analysis
//

include { DIFFHICAR } from '../../modules/local/bioc/diffhicar'

workflow DA {
    take:
    loops                   // peaks regions
    samplebedpe             // bedpe file per sample
    long_bedpe_postfix

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()

    ch_diffhicar = loops.map{meta, bin_size, peak -> [bin_size, peak]}
        .groupTuple()
        .cross(samplebedpe.map{[it[0].bin, it[1]]}.groupTuple())
        .map{ peak, long_bedpe ->
            [peak[0], peak[1].flatten(), long_bedpe[1].flatten()] }//bin_size, meta, peak, long_bedpe
        .groupTuple()
        .map{[it[0], it[1].flatten().unique(), it[2].flatten()]}
        .filter{it[1].size > 1} // filter by the bedpe files. Single bedpe means single group, no need to do differential analysis

    if(ch_diffhicar){
        switch(params.da_tool){
            case "edger":
                DIFFHICAR(ch_diffhicar, long_bedpe_postfix)
                ch_annotation_files = DIFFHICAR.out.anno
                ch_versions = ch_versions.mix(DIFFHICAR.out.versions.ifEmpty(null))
                ch_multiqc_files = ch_multiqc_files.mix(DIFFHICAR.out.stats.collect().ifEmpty(null))
                break
            default:
                DIFFHICAR(ch_diffhicar, long_bedpe_postfix)
                ch_annotation_files = DIFFHICAR.out.anno
                ch_versions = ch_versions.mix(DIFFHICAR.out.versions.ifEmpty(null))
                ch_multiqc_files = ch_multiqc_files.mix(DIFFHICAR.out.stats.collect().ifEmpty(null))
                break
        }
    }

    emit:
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
