//
// Call HiPeak
//

include { R1_PEAK } from './hipeak/fragment_peak'
include { HI_PEAK } from './hipeak/call_hipeak'

workflow HIPEAK {
    take:
    fasta                // channel: [ path(fasta) ]
    chrom_sizes          // channel: [ path(chrom_size) ]
    digest_genome        // channel: [ path(cut) ]
    mappability          // channel: [ path(mappability) ]
    gtf                  // channel: [ path(gtf) ]
    distalpair           // channel: [ val(meta), [pairs] ]
    hdf5                 // channel: [ val(meta), [hdf5] ]
    peak                 // channel: [ val(meta), [peak] ]
    r1_pval_thresh       // value
    short_bed_postfix    // value
    maps_3d_ext          // value

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_circos_files         = Channel.empty()
    ch_track_files          = Channel.empty()
    ch_annotation_files     = Channel.empty()

    R1_PEAK(
        distalpair,
        chrom_sizes,
        digest_genome,
        gtf,
        r1_pval_thresh,
        short_bed_postfix
    )
    ch_versions = ch_versions.mix(R1_PEAK.out.versions)
    //ch_trackfiles = ch_trackfiles.mix(
    //    R1_PEAK.out.bws.map{[it[0].id+"_R1",
    //        RelativePublishFolder.getPublishedFolder(workflow,
    //            'UCSC_BEDGRAPHTOBIGWIG_PER_R1_GROUP')+it[1].name]})

    // merge ATAC_PEAK with R1_PEAK by group id
    distalpair = hdf5.map{meta, bed -> [meta.group, bed]}
                                        .groupTuple()
    grouped_reads_peak = peak.map{[it[0].id, it[1]]}
                            .join(R1_PEAK.out.peak.map{[it[0].id, it[1]]})
                            .join(distalpair)
                            .map{[[id:it[0]], it[1], it[2], it[3]]}
    HI_PEAK(
        grouped_reads_peak,
        chrom_sizes,
        gtf,
        fasta,
        digest_genome,
        mappability,
        params.skip_peak_annotation,
        params.skip_diff_analysis,
        maps_3d_ext
    )
    ch_versions = ch_versions.mix(HI_PEAK.out.versions)
    //ch_trackfiles = ch_trackfiles.mix(
    //    HI_PEAK.out.bedpe
    //            .map{[it[0].id+"_HiPeak",
    //                RelativePublishFolder.getPublishedFolder(workflow,
    //                                    'ASSIGN_TYPE')+it[1].name]})

    ch_circos_files = HI_PEAK.out.bedpe

    emit:
    fragmentPeak    = R1_PEAK.out.peak
    bed4tfea        = HI_PEAK.out.bed4tfea
    circos          = ch_circos_files
    igv             = ch_track_files
    anno            = ch_annotation_files
    versions        = ch_versions          // channel: [ versions.yml ]
    mqc             = ch_multiqc_files
}
