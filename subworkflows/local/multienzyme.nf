/*
 * Createing Genomic Features Files
 */
include { initOptions } from './functions'
params.options = [:]
options        = initOptions(params.options)

include { MAPS_CUT                     } from '../maps/cut'                     addParams(options: options.maps_cut)
include { MAPS_FEND                    } from '../maps/fend'                    addParams(options: options.maps_fend)
include { GENMAP_MAPPABILITY           } from '../genmap/mappability'           addParams(options: options.genmap_mappability)
include { SEQLEVELS_STYLE              } from '../bioc/seqlevelsstyle'
include { ENSEMBL_UCSC_CONVERT
    ENSEMBL_UCSC_CONVERT as ENSEMBL_UCSC_CONVERT2       } from '../bioc/ensembl_ucsc_convert'        addParams(options: [args: "toUCSC", publish_dir:''])
include { UCSC_WIGTOBIGWIG             } from '../ucsc/wigtobigwig'             addParams(options: options.ucsc_wigtobigwig)
include { UCSC_BIGWIGAVERAGEOVERBED    } from '../ucsc/bigwigaverageoverbed'    addParams(options: options.maps_mapability)
include { MAPS_MERGE                   } from '../maps/merge'                   addParams(options: options.maps_merge)
include { MAPS_FEATURE                 } from '../maps/feature'                 addParams(options: options.maps_feature)

workflow MAPS_MULTIENZYME {
    take:
    fasta        // channel: [ path(fasta) ]
    cool_bin     // channel: [ val(bin) ]
    chromsizes   // channel: [ path(chromsizes) ]

    main:
    ch_version = MAPS_CUT(fasta, cool_bin).version
    MAPS_FEND(MAPS_CUT.out.cut, chromsizes)
    if(!params.mappability){
        GENMAP_MAPPABILITY(fasta)
        ch_version = ch_version.mix(GENMAP_MAPPABILITY.out.version)
        mappability = UCSC_WIGTOBIGWIG(GENMAP_MAPPABILITY.out.wig, chromsizes).bw
        ch_version = ch_version.mix(UCSC_WIGTOBIGWIG.out.version)
    }else{
        mappability = Channel.fromPath(params.mappability, checkIfExists: true)
    }
    seqlevelsstyle = SEQLEVELS_STYLE(MAPS_FEND.out.bed.map{it[1]}.collect().map{it[0]}).seqlevels_style
    if("$seqlevelsstyle" != "UCSC"){
        ENSEMBL_UCSC_CONVERT(MAPS_FEND.out.bed)
        ENSEMBL_UCSC_CONVERT2(mappability.map{[cool_bin, it]})
        ch_version = ch_version.mix(ENSEMBL_UCSC_CONVERT.out.version)
        UCSC_BIGWIGAVERAGEOVERBED(ENSEMBL_UCSC_CONVERT.out.tab, ENSEMBL_UCSC_CONVERT2.out.tab.map{it[1]})
    }else{
        UCSC_BIGWIGAVERAGEOVERBED(MAPS_FEND.out.bed, mappability)
    }
    ch_version = ch_version.mix(UCSC_BIGWIGAVERAGEOVERBED.out.version)
    MAPS_MERGE(MAPS_CUT.out.cut.join(UCSC_BIGWIGAVERAGEOVERBED.out.tab))

    MAPS_FEATURE(MAPS_MERGE.out.map, chromsizes)

    emit:
    bin_feature              = MAPS_FEATURE.out.bin_feature      // channel: [ val(bin_size), path(bin_feature) ]
    version                  = ch_version                        // channel: [ path(version) ]
}
