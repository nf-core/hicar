//
// differential analysis
//

include { DIFFHICAR } from '../../modules/local/bioc/diffhicar'
include { DIFFHIC   } from '../../modules/local/bioc/diffhic'
include { DIFFSET   } from '../../modules/local/bioc/diffset'
include { HICEXPLORER_DIFFHIC } from './diff_ana/chicdiff'

workflow DA {
    take:
    loops                   // peaks regions
    samplebedpe             // bedpe file per sample

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty() // TODO
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_diff_files           = Channel.empty()

    if(params.da_tool!='hicexplorer'){
        ch_diffhicar = loops.map{meta, bin_size, peak -> [bin_size, peak]}
            .groupTuple()
            .cross(samplebedpe.map{[it[0].bin, it[1]]}.groupTuple())
            .map{ peak, long_bedpe ->
                [peak[0], peak[1].flatten(), long_bedpe[1].flatten()] }//bin_size, meta, peak, long_bedpe
            .groupTuple()
            .map{[it[0], it[1].flatten().unique(), it[2].flatten()]}
            .filter{it[1].size > 1} // filter by the bedpe files. Single bedpe means single group, no need to do differential analysis
    }else{
        ch_diffhicar = loops
    }

    if(ch_diffhicar){
        switch(params.da_tool){
            case "edger":
                DIFFHICAR(ch_diffhicar, params.long_bedpe_postfix)
                ch_annotation_files = DIFFHICAR.out.anno
                ch_versions = ch_versions.mix(DIFFHICAR.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(DIFFHICAR.out.stats.collect().ifEmpty(null))
                ch_diff_files = ch_diff_files.mix(DIFFHICAR.out.diff)
                break
            case "diffhic":
                DIFFHIC(ch_diffhicar, params.long_bedpe_postfix)
                ch_annotation_files = DIFFHIC.out.anno
                ch_versions = ch_versions.mix(DIFFHIC.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(DIFFHIC.out.stats.collect().ifEmpty(null))
                ch_diff_files = ch_diff_files.mix(DIFFHIC.out.diff)
                break
            case "setOperation":
                DIFFSET(ch_diffhicar, params.long_bedpe_postfix)
                ch_annotation_files = DIFFSET.out.anno
                ch_versions = ch_versions.mix(DIFFSET.out.versions)
                ch_diff_files = ch_diff_files.mix(DIFFSET.out.diff)
                break
            case "hicexplorer":
                HICEXPLORER_DIFFHIC(ch_diffhicar)
                ch_versions = ch_versions.mix(HICEXPLORER_DIFFHIC.out.versions)
                ch_diff_files = ch_diff_files.mix(HICEXPLORER_DIFFHIC.out.diff)
                break
            default:
                DIFFHIC(ch_diffhicar, params.long_bedpe_postfix)
                ch_annotation_files = DIFFHIC.out.anno
                ch_versions = ch_versions.mix(DIFFHIC.out.versions)
                ch_multiqc_files = ch_multiqc_files.mix(DIFFHIC.out.stats.collect().ifEmpty(null))
                break
        }
    }

    emit:
    diff            = ch_diff_files
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
