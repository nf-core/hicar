/*
 * call APA by cooltools
 */

include { COOLTOOLS_PILEUP                 } from '../../../modules/local/cooltools/pileup'
include { PLOTNPZ_BY_COOLTOOLS                } from '../../../modules/local/cooltools/plotnpz'

workflow COOLTOOLS_APACALLER {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    peaks                  // path(1D peaks)
    format                 // output plot format

    main:
    COOLTOOLS_PILEUP(matrix, peaks)
    ch_version = COOLTOOLS_PILEUP.out.versions

    PLOTNPZ_BY_COOLTOOLS(COOLTOOLS_PILEUP.out.npz, format)


    emit:
    plot      = PLOTNPZ_BY_COOLTOOLS.out.plot                    // channel: [ val(meta), path(plot)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
