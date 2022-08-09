/*
 * Call interaction peaks by MAPS
 */

include { MAPS_MULTIENZYME } from './maps_multienzyme'
include { MAPS_PEAK        } from './maps_peak'

workflow MAPS {
    take:
    fasta                           // [ genome fa ]
    chrom_sizes                     // [ chromsizes ]
    mappability                     // [ bigwig file ]
    bin_size                        // values: bin size, 5000, 10000
    site                            // values: eg. 'GATC 1'
    reads                           // [ meta, [bedgraph] ] cooler dump reads for each group
    mergedpeak                      // [ peaks ] merged bed file
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome reads for group
    merge_map_py_source             // scripts
    feature_frag2bin_source         // scripts
    make_maps_runfile_source        // scripts
    juicer_tools                    // scripts
    long_bedpe_postfix              // values
    short_bed_postfix               // values
    maps_3d_ext                     // values

    main:
    // generate features
    background = MAPS_MULTIENZYME(fasta,
                    bin_size,
                    chrom_sizes,
                    mappability,
                    merge_map_py_source,
                    feature_frag2bin_source,
                    params.enzyme,
                    site).bin_feature
    ch_versions = MAPS_MULTIENZYME.out.versions

    reads_peak   = reads
                    .map{ meta, reads ->
                            [meta.id, reads]} // here id is group
                    .combine(mergedpeak)// group, reads, peaks
                    .cross(bedpe.map{[it[0].id, it[0].bin, it[1]]})// group, bin, bedpe
                    .map{ short_bed, long_bedpe -> //[bin_size, group, macs2, long_bedpe, short_bed]
                            [long_bedpe[1], short_bed[0], short_bed[2], long_bedpe[2], short_bed[1]]}
    background.cross(reads_peak)
                .map{ background, reads -> //[group, bin_size, macs2, long_bedpe, short_bed, background]
                        [[id:reads[1]], background[0], reads[2], reads[3], reads[4], background[1]]}
                .set{ maps_input }
    MAPS_PEAK(
        maps_input,
        make_maps_runfile_source,
        chrom_sizes,
        params.juicer_jvm_params,
        juicer_tools,
        long_bedpe_postfix,
        short_bed_postfix,
        maps_3d_ext)
    ch_versions = ch_versions.mix(MAPS_PEAK.out.versions.ifEmpty(null))
    ch_multiqc_files = MAPS_PEAK.out.stats.collect().ifEmpty(null)

    emit:
    interactions = MAPS_PEAK.out.peak           // channel: [ meta, bin_size, path(bedpe) ]
    mqc          = ch_multiqc_files             // channel: [ path(mqc) ]
    versions     = ch_versions                  // channel: [ path(version) ]
}
