/*
 * Call interaction peaks by MAPS
 */

include { CIRCOS_PREPARE            } from '../../modules/local/circos/circos_prepare'
include { CIRCOS                    } from '../../modules/local/circos/circos'

workflow RUN_CIRCOS {
    take:
    bedpe            // channel: [ path(bedpe) ]
    gtf              // channel: [ path(gtf) ]
    chromsize        // channel: [ path(chromsize) ]
    ucscname         // channel: [ val(ucscname) ]
    config           // channel: [ path(config) ]

    main:
    //create circos config
    //input=path(bedpe), val(ucscname), path(gtf), path(chromsize)
    ch_version = CIRCOS_PREPARE(bedpe.combine(ucscname).combine(gtf).combine(chromsize)).versions.first()
    //plot
    CIRCOS(CIRCOS_PREPARE.out.circos, config)
    ch_version = ch_version.mix(CIRCOS.out.versions)

    emit:
    circos       = CIRCOS.out.circos            // channel: [ path(png) ]
    versions     = ch_version                   // channel: [ path(version) ]
}
