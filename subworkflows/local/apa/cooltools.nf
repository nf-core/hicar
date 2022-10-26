/*
 * call APA by cooltools
 */

include { COOLTOOLS_PILEUP                 } from '../../../modules/local/cooltools/pileup'
include { COOLTOOLS_PLOTNPZ                } from '../../../modules/local/cooltools/plotnpz'

workflow COOLTOOLS_APACALLER {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    peaks                  // path(1D peaks)

    main:
    COOLTOOLS_PILEUP(matrix, peaks).npz | COOLTOOLS_PLOTNPZ

    ch_version = COOLTOOLS_PILEUP.out.versions.ifEmpty(null)

    emit:
    png       = COOLTOOLS_PLOTNPZ.out.png                        // channel: [ val(meta), path(pngs)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
