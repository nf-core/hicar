/*
 * Call interaction peaks by MAPS
 */

include { MAPS_MULTIENZYME } from './maps_multienzyme'
include { MAPS_PEAK        } from './maps_peak'

workflow MAPS {
    take:
    bedpe                           // 2D: [ val(meta), [bedpe] ] merged intra_chromosome reads for group
    reads                           // 1D: [ meta, [group_reads_count_bedgraph], [merged_peaks_bed] ]
    additional_param                // [
                                    //    0 values: cuting_site, output of prepare genome
                                    //    1 [ genome fa ], output of prepare genome
                                    //    2 [ chrom_size ], output of prepare genome
                                    //    3 [ mappability bigwig ], output of prepare genome
                                    // ]

    main:
    // generate features
    background = MAPS_MULTIENZYME(
        bedpe.map{ [ it[0].bin ] }.unique(),                   // bin_size
        additional_param.map{[ it[0], it[1], it[2], it[3] ]}   // site, fasta, chrom_size, mappability_bw
        ).bin_feature
    ch_versions = MAPS_MULTIENZYME.out.versions

    reads_peak   = reads
                    .map{ meta, reads, peaks ->
                            [meta.id, reads, peaks]} // here id is group
                    .cross(bedpe.map{[it[0].id, it[0].bin, it[1]]})// group, bin, bedpe
                    .map{ short_bed, long_bedpe -> //[bin_size, group, macs2, long_bedpe, short_bed]
                            [long_bedpe[1], short_bed[0], short_bed[2], long_bedpe[2], short_bed[1]]}
    background.cross(reads_peak)
                .map{ background, reads -> //[group, bin_size, macs2, long_bedpe, short_bed, background]
                        [[id:reads[1]], background[0], reads[2], reads[3], reads[4], background[1]]}
                .set{ maps_input }
    MAPS_PEAK(
        maps_input,
        additional_param.map{[ it[1] ]},            //chrom_size
        )
    ch_versions = ch_versions.mix(MAPS_PEAK.out.versions)
    ch_multiqc_files = MAPS_PEAK.out.stats.collect().ifEmpty(null)

    emit:
    interactions = MAPS_PEAK.out.peak           // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
