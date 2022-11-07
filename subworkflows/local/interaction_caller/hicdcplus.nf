/*
 * Call interaction peaks by HiC-CD+
 */

include { HICDCPLUS_FEATURES             } from '../../../modules/local/hicdcplus/features'
include { HICDCPLUS_CALLLOOPS            } from '../../../modules/local/hicdcplus/callloops'

workflow HICDCPLUS {
    take:
    fasta                           // [ genome fa ]
    chrom_sizes                     // [ chromsizes ]
    mappability                     // [ bigwig file ]
    bin_size                        // values: bin size, 5000, 10000
    site                            // values: eg. 'GATC 1'
    reads                           // [ meta, [bedgraph] ] cooler dump reads for each group
    mergedpeak                      // [ peaks ] merged bed file
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome reads for group

    main:
    ch_multiqc_files = Channel.empty()
    ch_loop = Channel.empty()

    // generate features
    ch_versions = HICDCPLUS_FEATURES(
                    bin_size.combine(site)
                        .combine(fasta)
                        .combine(chrom_sizes)
                        .combine(mappability)
                    ).versions
    // call loops
    ch_reads = reads.combine(HICDCPLUS_FEATURES.out.features).map{[[id:it[0].id, bin:it[2]], it[1], it[3]]}
        .combine(bedpe, by: 0)
    ch_reads.view()
    ch_loop = HICDCPLUS_CALLLOOPS(ch_reads, chrom_sizes).interactions
    ch_versions = ch_versions.mix(HICDCPLUS_CALLLOOPS.out.versions.ifEmpty([]))

    emit:
    interactions = ch_loop                      // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
