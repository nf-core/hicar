/*
 * Call interaction peaks by MAPS
 */
params.options = [:]

include { MAPS_MAPS             } from '../../modules/local/maps/maps'         addParams(options: params.options.maps_maps)
include { MAPS_CALLPEAK         } from '../../modules/local/maps/callpeak'     addParams(options: params.options.maps_callpeak)
include { READS_SUMMARY
    as MAPS_STATS               } from '../../modules/local/reads_summary'      addParams(options: params.options.maps_stats)
include { MAPS_REFORMAT         } from '../../modules/local/maps/reformat'     addParams(options: params.options.maps_reformat)

workflow MAPS_PEAK {
    take:
    reads        // channel: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]

    main:
    //create parameter table
    //input=val(meta), val(bin_size), path(macs2), path(long_bedpe), path(short_bed), path(background)
    //maps from bedpe
    ch_version = MAPS_MAPS(reads).versions
    //regression and peak calling
    peak = MAPS_CALLPEAK(MAPS_MAPS.out.maps).peak
    ch_version = ch_version.mix(MAPS_CALLPEAK.out.versions)
    //peak formatting
    MAPS_REFORMAT(peak)
    ch_version = ch_version.mix(MAPS_REFORMAT.out.versions)
    //merge stats
    MAPS_STATS(MAPS_CALLPEAK.out.summary.map{it[2]}.collect())

    emit:
    peak         = MAPS_REFORMAT.out.bedpe      // channel: [ path(bedpe) ]
    stats        = MAPS_STATS.out.summary       // channel: [ path(stats) ]
    versions     = ch_version                   // channel: [ path(version) ]
}
