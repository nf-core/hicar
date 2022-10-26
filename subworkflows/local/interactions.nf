//
// Call interaction loops
// Ask all the output of each module must contain channel:
// versions: [path]; interactions: [ meta, bin_size, [bedpe] ]
// optional output: mqc: [path];
//

include { HICDCPLUS          } from './interaction_caller/hicdcplus'
include { MAPS               } from './interaction_caller/maps'
include { MERGE_INTERACTIONS } from '../../modules/local/bioc/merge_interactions'

workflow INTERACTIONS {
    take:
    // 1D peaks
    reads                           // [ meta, [bedgraph] ] cooler dump ATAC reads for each group
    mergedpeak                      // [ peaks ] merged bed file
    // 2D signals
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome and inter_chromosome reads for group
    tagdir                          // [ val(meta), [tagdir] ]
    // reference
    fasta                           // [ genome fa ]
    chrom_sizes                     // [ chromsizes ]
    mappability                     // [ bigwig file ]
    site                            // values: eg. 'GATC 1'
    genome                          // value: UCSC genome name: hg38, mm10 ...
    // resolutions
    bin_size                        // values: bin size, 5000, 10000
    // source
    merge_map_py_source             // scripts
    feature_frag2bin_source         // scripts
    make_maps_runfile_source        // scripts
    // other constant values
    long_bedpe_postfix              // values
    short_bed_postfix               // values
    maps_3d_ext                     // values

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()
    ch_loops                = Channel.empty() // a bed files channel

    switch(params.interactions_tool){
        case "maps":
            MAPS(
                fasta,
                chrom_sizes,
                mappability,
                bin_size,
                site,
                reads,
                mergedpeak,
                bedpe,
                merge_map_py_source,
                feature_frag2bin_source,
                make_maps_runfile_source,
                long_bedpe_postfix,
                short_bed_postfix,
                maps_3d_ext
            )
            ch_loops = MAPS.out.interactions
            ch_versions = MAPS.out.versions
            ch_annotation_files = MAPS.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = MAPS.out.interactions.map{
                meta, bin_size, bedpe ->
                    meta.bin = bin_size
                    [meta, bedpe]
            }
            break
        case "hicdcplus":
            HICDCPLUS(
                fasta,
                chrom_sizes,
                mappability,
                bin_size,
                site,
                reads,
                mergedpeak,
                bedpe
            )
            ch_loops = HICDCPLUS.out.interactions
            ch_versions = HICDCPLUS.out.versions
            ch_annotation_files = HICDCPLUS.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = HICDCPLUS.out.interactions.map{
                meta, bin_size, bedpe ->
                    meta.bin = bin_size
                    [meta, bedpe]
            }
            break
        default:
            MAPS(
                fasta,
                chrom_sizes,
                mappability,
                bin_size,
                site,
                reads,
                mergedpeak,
                bedpe,
                merge_map_py_source,
                feature_frag2bin_source,
                make_maps_runfile_source,
                juicer_tools,
                long_bedpe_postfix,
                short_bed_postfix,
                maps_3d_ext
            )
            ch_loops = MAPS.out.interactions
            ch_versions = MAPS.out.versions
            ch_annotation_files = MAPS.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = MAPS.out.interactions.map{
                meta, bin_size, bedpe ->
                    meta.bin = bin_size
                    [meta, bedpe]
            }
            break
    }
    // merge loops
    MERGE_INTERACTIONS(
        ch_loops.map{
            meta, bin_size, interactions -> [bin_size, interactions]
        }.groupTuple()
    )
    ch_versions = ch_versions.mix(MERGE_INTERACTIONS.out.versions.ifEmpty(null))

    emit:
    loops           = ch_loops                              // channel: [ meta, bin_size, path(bedpe) ]
    mergedloops     = MERGE_INTERACTIONS.out.interactions   // channel: [ bin_siz, path(bedpe) ]
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions                           // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
