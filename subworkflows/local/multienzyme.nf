/*
 * Createing Genomic Features Files
 */
params.options = [:]

include { MAPS_CUT                     } from '../../modules/local/maps/cut'                     addParams(options: params.options.maps_cut)
include { MAPS_FEND                    } from '../../modules/local/maps/fend'                    addParams(options: params.options.maps_fend)
include { GENMAP_INDEX                 } from '../../modules/nf-core/modules/genmap/index/main'           addParams(options: params.options.genmap_index)
include { GENMAP_MAPPABILITY           } from '../../modules/nf-core/modules/genmap/mappability/main'           addParams(options: params.options.genmap_mappability)
include { SEQLEVELS_STYLE              } from '../../modules/local/bioc/seqlevelsstyle'
include { ENSEMBL_UCSC_CONVERT
    ENSEMBL_UCSC_CONVERT as ENSEMBL_UCSC_CONVERT2       } from '../../modules/local/bioc/ensembl_ucsc_convert'        addParams(options: params.options.ensembl_ucsc_convert)
include { UCSC_WIGTOBIGWIG             } from '../../modules/nf-core/modules/ucsc/wigtobigwig/main'             addParams(options: params.options.ucsc_wigtobigwig)
include { UCSC_BIGWIGAVERAGEOVERBED    } from '../../modules/nf-core/modules/ucsc/bigwigaverageoverbed/main'    addParams(options: params.options.maps_mapability)
include { MAPS_MERGE                   } from '../../modules/local/maps/merge'                   addParams(options: params.options.maps_merge)
include { MAPS_FEATURE                 } from '../../modules/local/maps/feature'                 addParams(options: params.options.maps_feature)

workflow MAPS_MULTIENZYME {
    take:
    fasta        // channel: [ path(fasta) ]
    cool_bin     // channel: [ val(bin) ]
    chromsizes   // channel: [ path(chromsizes) ]

    main:
    if(params.maps_digest_file && params.enzyme.toLowerCase() != "mnase"){
        ch_version = Channel.empty()
        ch_digest = cool_bin.combine(Channel.fromPath(params.maps_digest_file))
        MAPS_FEND(ch_digest, chromsizes)
    }else{
        ch_version = MAPS_CUT(fasta, cool_bin).version
        ch_digest = MAPS_CUT.out.cut
        MAPS_FEND(ch_digest, chromsizes)
    }
    if(!params.mappability){
        GENMAP_INDEX(fasta).index | GENMAP_MAPPABILITY
        ch_version = ch_version.mix(GENMAP_MAPPABILITY.out.version)
        mappability = UCSC_WIGTOBIGWIG(GENMAP_MAPPABILITY.out.wig, chromsizes).bw
        ch_version = ch_version.mix(UCSC_WIGTOBIGWIG.out.version)
    }else{
        mappability = Channel.fromPath(params.mappability, checkIfExists: true)
    }
    seqlevelsstyle = SEQLEVELS_STYLE(MAPS_FEND.out.bed.map{it[1]}.collect().map{it[0]}).seqlevels_style
    if("$seqlevelsstyle" != "UCSC"){
        ENSEMBL_UCSC_CONVERT(MAPS_FEND.out.bed)
        ENSEMBL_UCSC_CONVERT2(cool_bin.combine(mappability))
        ch_version = ch_version.mix(ENSEMBL_UCSC_CONVERT.out.version)
        UCSC_BIGWIGAVERAGEOVERBED(ENSEMBL_UCSC_CONVERT.out.tab.map{[['id':'background', 'bin_size':it[0]], it[1]]},
                                    ENSEMBL_UCSC_CONVERT2.out.tab.map{it[1]})
    }else{
        UCSC_BIGWIGAVERAGEOVERBED(MAPS_FEND.out.bed.map{[['id':'background', 'bin_size':it[0]], it[1]]}, mappability)
    }
    ch_version = ch_version.mix(UCSC_BIGWIGAVERAGEOVERBED.out.version)
    MAPS_MERGE(ch_digest.cross(UCSC_BIGWIGAVERAGEOVERBED.out.tab.map{[it[0].bin_size, it[1]]}).map{[it[0][0], it[0][1], it[1][1]]})

    MAPS_FEATURE(MAPS_MERGE.out.map, chromsizes)

    emit:
    mappability              = mappability                       // channel: [ path(bw) ]
    bin_feature              = MAPS_FEATURE.out.bin_feature      // channel: [ val(bin_size), path(bin_feature) ]
    version                  = ch_version                        // channel: [ path(version) ]
}
