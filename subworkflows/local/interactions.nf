//
// Call interaction loops
//

include { HICDCPLUS       } from './interaction_caller/hicdcplus'
include { MAPS            } from './interaction_caller/maps'
include { HOMER           } from './interaction_caller/homer'

workflow INTERACTIONS {
    take:
    fasta                           // [ genome fa ]
    chrom_sizes                     // [ chromsizes ]
    mappability                     // [ bigwig file ]
    bin_size                        // values: bin size, 5000, 10000
    site                            // values: eg. 'GATC 1'
    reads                           // [ meta, [bedgraph] ] cooler dump ATAC reads for each group
    mergedpeak                      // [ peaks ] merged bed file
    bedpe                           // [ val(meta), [bedpe] ] merged intra_chromosome and tnter_chromosome reads for group
    tagdir                          // [ val(meta), [tagdir] ]
    merge_map_py_source             // scripts
    feature_frag2bin_source         // scripts
    make_maps_runfile_source        // scripts
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
        case "homer":
            HOMER(
                fasta,
                chrom_sizes,
                mappability,
                bin_size,
                site,
                tagdir,
                mergedpeak,
                bedpe
            )
            ch_loops = HOMER.out.interactions
            ch_versions = HOMER.out.versions
            ch_annotation_files = HOMER.out.interactions.map{
                meta, bin_size, interactions -> [meta.id+bin_size, interactions]}
            ch_circos_files = HOMER.out.interactions.map{
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


    emit:
    loops           = ch_loops             // channel: [ meta, bin_size, path(bedpe) ]
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
