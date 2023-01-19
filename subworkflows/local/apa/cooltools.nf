/*
 * call APA by cooltools
 */

include { COOLTOOLS_PILEUP                 } from '../../../modules/local/cooltools/pileup'
include { PLOTNPZ_BY_COOLTOOLS                } from '../../../modules/local/cooltools/plotnpz'

workflow COOLTOOLS_APACALLER {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    peaks                  // path(1D peaks)

    main:
    COOLTOOLS_PILEUP(matrix, peaks).npz | PLOTNPZ_BY_COOLTOOLS

    ch_version = COOLTOOLS_PILEUP.out.versions.ifEmpty(null)

    emit:
    png       = PLOTNPZ_BY_COOLTOOLS.out.png                        // channel: [ val(meta), path(pngs)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
