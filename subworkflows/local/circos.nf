/*
 * Call interaction peaks by MAPS
 */
params.options = [:]

include { PREPARE_CIRCOS            } from '../../modules/local/circos/prepare_circos'    addParams(options: params.options.prepare_circos)
include { CIRCOS                    } from '../../modules/local/circos/circos'            addParams(options: params.options.circos)

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
    ch_version = PREPARE_CIRCOS(bedpe.combine(ucscname).combine(gtf).combine(chromsize)).version
    //plot
    CIRCOS(PREPARE_CIRCOS.out.circos.combine(config))
    ch_version = ch_version.mix(CIRCOS.out.version)

    emit:
    circos       = CIRCOS.out.circos            // channel: [ path(png) ]
    version      = ch_version                   // channel: [ path(version) ]
}
