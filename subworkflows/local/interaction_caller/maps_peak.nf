/*
 * Call interaction peaks by MAPS
 */

include { MAPS_MAPS             } from '../../../modules/local/maps/maps'
include { MAPS_CALLPEAK         } from '../../../modules/local/maps/callpeak'
include { READS_SUMMARY
    as MAPS_STATS               } from '../../../modules/local/reads_summary'
include { MAPS_REFORMAT         } from '../../../modules/local/maps/reformat'
include { MAPS_RAW2BG2          } from '../../../modules/local/maps/raw2bg2'
include { COOLER_LOAD           } from '../../../modules/local/cooler/load/main'
include { COOLER_MERGE          } from '../../../modules/nf-core/cooler/merge/main'
include { COOLER_ZOOMIFY
    as COOLER_ZOOMIFY_MAPS      } from '../../../modules/nf-core/cooler/zoomify/main'

workflow MAPS_PEAK {
    take:
    reads                     // channel: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]
    chromsizes                // channel: [ path(chromsizes) ]

    main:
    //create parameter table
    //input=val(meta), val(bin_size), path(macs2), path(long_bedpe), path(short_bed), path(background)
    //maps from bedpe
    ch_version = MAPS_MAPS(reads, params.long_bedpe_postfix, params.short_bed_postfix).versions
    //regression and peak calling
    peak = MAPS_CALLPEAK(MAPS_MAPS.out.maps).peak
    ch_version = ch_version.mix(MAPS_CALLPEAK.out.versions)
    if(params.create_maps_signal){
        //create cooler files for raw read matrix
        MAPS_RAW2BG2(MAPS_CALLPEAK.out.signal)
        ch_version = ch_version.mix(MAPS_RAW2BG2.out.versions)
        COOLER_LOAD(MAPS_RAW2BG2.out.bg2, chromsizes)
        ch_version = ch_version.mix(COOLER_LOAD.out.versions)
        // Merge contacts
        COOLER_LOAD.out.cool
                    .map{
                        meta, bin, cool ->
                        [[id:meta.id, bin:bin], cool]
                    }
                    .set{ch_cooler}
        COOLER_MERGE(ch_cooler)
        // create mcooler file for visualization
        COOLER_ZOOMIFY_MAPS(COOLER_MERGE.out.cool)
    }

    //peak formatting
    MAPS_REFORMAT(peak, params.maps_3d_ext)
    ch_version = ch_version.mix(MAPS_REFORMAT.out.versions)
    //merge stats
    MAPS_STATS(MAPS_CALLPEAK.out.summary.map{it[2]}.collect())

    emit:
    peak         = MAPS_REFORMAT.out.bedpe      // channel: [ meta, bin_size, path(bedpe) ]
    stats        = MAPS_STATS.out.summary       // channel: [ path(stats) ]
    versions     = ch_version                   // channel: [ path(version) ]
}
