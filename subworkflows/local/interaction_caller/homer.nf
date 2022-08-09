/*
 * Call interaction peaks by HiC-CD+
 */

include { HOMER_ANALYZEHIC             } from '../../../modules/local/homer/analyze_hic'

workflow HOMER {
    take:
    fasta                           // [ genome fa ]
    chrom_sizes                     // [ chromsizes ]
    mappability                     // [ bigwig file ]
    bin_size                        // values: bin size, 5000, 10000
    site                            // values: eg. 'GATC 1'
    tagdir                          // [ meta, [tagdir] ] cooler dump ATAC reads for each group
    mergedpeak                      // [ peaks ] merged bed file
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome and tnter_chromosome reads for group

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()
    ch_versions = Channel.empty()

    // call interactions
    tagdir.view()
    tag = tagdir.combine(bin_size).map{
        meta, tag, bin -> [ meta, bin, tag]
    }
    ch_loop = HOMER_ANALYZEHIC(tag).interactions
    ch_versions = HOMER_ANALYZEHIC.out.versions


    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
