/*
 * Createing Genomic Features Files
 */

include { BIOC_ENZYMECUT               } from '../../modules/local/bioc/enzyme_cut'
include { MAPS_CUT                     } from '../../modules/local/maps/cut'
include { MAPS_FEND                    } from '../../modules/local/maps/fend'
include { GENMAP_INDEX                 } from '../../modules/nf-core/modules/genmap/index/main'
include { GENMAP_MAPPABILITY           } from '../../modules/nf-core/modules/genmap/mappability/main'
include { SEQLEVELS_STYLE              } from '../../modules/local/bioc/seqlevelsstyle'
include { ENSEMBL_UCSC_CONVERT
    ENSEMBL_UCSC_CONVERT as ENSEMBL_UCSC_CONVERT2       } from '../../modules/local/bioc/ensembl_ucsc_convert'
include { UCSC_WIGTOBIGWIG             } from '../../modules/nf-core/modules/ucsc/wigtobigwig/main'
include { UCSC_BIGWIGAVERAGEOVERBED    } from '../../modules/nf-core/modules/ucsc/bigwigaverageoverbed/main'
include { MAPS_MERGE                   } from '../../modules/local/maps/merge'
include { MAPS_FEATURE                 } from '../../modules/local/maps/feature'

workflow MAPS_MULTIENZYME {
    take:
    fasta                      // channel: [ path(fasta) ]
    cool_bin                   // channel: [ val(bin) ]
    chromsizes                 // channel: [ path(chromsizes) ]
    merge_map_py_source        // channel: [ file(merge_map_py_source) ]
    feature_frag2bin_source    // channel: [ file(feature_frag2bin_source) ]
    enzyme                     // values
    site                       // values

    main:
    if(params.maps_digest_file && params.enzyme.toLowerCase() != "mnase"){
        ch_version = Channel.empty()
        ch_digest = cool_bin.combine(Channel.fromPath(params.maps_digest_file))
        MAPS_FEND(ch_digest, chromsizes)
    }else{
        if(params.enzyme.toLowerCase() != "mnase"){
            ch_version = BIOC_ENZYMECUT(fasta, cool_bin, site, enzyme).versions
            ch_digest = BIOC_ENZYMECUT.out.cut
        }else{
            ch_version = MAPS_CUT(fasta, cool_bin, site, enzyme).versions
            ch_digest = MAPS_CUT.out.cut
        }
        MAPS_FEND(ch_digest, chromsizes)
    }
    if(!params.mappability){
        GENMAP_INDEX(fasta).index | GENMAP_MAPPABILITY
        ch_version = ch_version.mix(GENMAP_MAPPABILITY.out.versions)
        mappability = UCSC_WIGTOBIGWIG(GENMAP_MAPPABILITY.out.wig.map{[[id:'mappability'], it]}, chromsizes).bw.map{it[1]}
        ch_version = ch_version.mix(UCSC_WIGTOBIGWIG.out.versions)
    }else{
        mappability = Channel.fromPath(params.mappability, checkIfExists: true)
    }
    seqlevelsstyle = SEQLEVELS_STYLE(MAPS_FEND.out.bed.map{it[1]}.collect().map{it[0]}).seqlevels_style
    if("$seqlevelsstyle" != "UCSC"){
        ENSEMBL_UCSC_CONVERT(MAPS_FEND.out.bed)
        ENSEMBL_UCSC_CONVERT2(cool_bin.combine(mappability))
        ch_version = ch_version.mix(ENSEMBL_UCSC_CONVERT.out.versions)
        UCSC_BIGWIGAVERAGEOVERBED(ENSEMBL_UCSC_CONVERT.out.tab.map{[['id':'background', 'bin_size':it[0]], it[1]]},
                                    ENSEMBL_UCSC_CONVERT2.out.tab.map{it[1]})
    }else{
        UCSC_BIGWIGAVERAGEOVERBED(MAPS_FEND.out.bed.map{[['id':'background', 'bin_size':it[0]], it[1]]}, mappability)
    }
    ch_version = ch_version.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions)
    MAPS_MERGE(ch_digest.cross(UCSC_BIGWIGAVERAGEOVERBED.out.tab.map{[it[0].bin_size, it[1]]}).map{[it[0][0], it[0][1], it[1][1]]},
                merge_map_py_source)

    MAPS_FEATURE(MAPS_MERGE.out.map, chromsizes, feature_frag2bin_source)

    emit:
    mappability              = mappability                       // channel: [ path(bw) ]
    bin_feature              = MAPS_FEATURE.out.bin_feature      // channel: [ val(bin_size), path(bin_feature) ]
    versions                 = ch_version                        // channel: [ path(version) ]
}
