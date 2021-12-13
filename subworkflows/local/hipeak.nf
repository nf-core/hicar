/*
 * Call interaction peaks by MAPS
 */
params.options = [:]

include { PREPARE_COUNTS            } from '../../modules/local/hipeak/prepare_counts'    addParams(options: params.options.parepare_counts)
include { CALL_HIPEAK               } from '../../modules/local/hipeak/callpeak'          addParams(options: params.options.call_hipeak)
include { ASSIGN_TYPE               } from '../../modules/local/hipeak/assign_type'       addParams(options: params.options.assign_type)
include { DIFF_HIPEAK               } from '../../modules/local/hipeak/diff_hipeak'       addParams(options: params.options.diff_hipeak)
include { BIOC_CHIPPEAKANNO
    as BIOC_CHIPPEAKANNO_HIPEAK     } from '../../modules/local/bioc/chippeakanno'        addParams(options: params.options.chippeakanno_hipeak)
include { BIOC_CHIPPEAKANNO
    as BIOC_CHIPPEAKANNO_DIFFHIPEAK } from '../../modules/local/bioc/chippeakanno'        addParams(options: params.options.chippeakanno_diffhipeak)
include { PAIR2BAM                  } from '../../modules/local/bioc/pair2bam'            addParams(options: params.options.pair2bam)

workflow HI_PEAK {
    take:
    peaks                // channel: [ meta, r2peak, r1peak, distalpair ]
    gtf                  // channel: [ path(gtf) ]
    fasta                // channel: [ path(fasta) ]
    digest_genome        // channel: [ path(digest_genome) ]
    mappability          // channel: [ path(mappability) ]
    skip_peak_annotation // value: params.skip_peak_annotation
    skip_diff_analysis   // value: params.skip_diff_analysis


    main:
    //create count table
    //input=val(meta), path(r2peak), path(r1peak), path(distalpair)
    ch_version = PREPARE_COUNTS(peaks.combine(mappability).combine(fasta).combine(digest_genome)).versions
    //regression and peak calling
    CALL_HIPEAK(PREPARE_COUNTS.out.counts)
    ch_version = ch_version.mix(CALL_HIPEAK.out.versions)
    PAIR2BAM(CALL_HIPEAK.out.peak.join(peaks.map{[it[0], it[3]]}))

    //assign type for peak
    ASSIGN_TYPE(CALL_HIPEAK.out.peak)
    ch_version = ch_version.mix(ASSIGN_TYPE.out.versions)
    // annotation
    if(!skip_peak_annotation){
        BIOC_CHIPPEAKANNO_HIPEAK(ASSIGN_TYPE.out.peak
                                            .map{it[1]}
                                            .filter{ it.readLines().size > 1 }
                                            .collect()
                                            .map{["HiPeak", it]}, gtf)
        ch_version = ch_version.mix(BIOC_CHIPPEAKANNO_HIPEAK.out.versions.ifEmpty(null))
    }
    //differential analysis
    stats = Channel.empty()
    diff = Channel.empty()
    if(!skip_diff_analysis){
        hipeaks = ASSIGN_TYPE.out.peak.map{it[1]}.collect().filter{it.size()>1}
        if(hipeaks){
            DIFF_HIPEAK(hipeaks,
                        peaks.map{it[3]}.collect())
            ch_version = ch_version.mix(DIFF_HIPEAK.out.versions.ifEmpty(null))
            stats = DIFF_HIPEAK.out.stats
            diff = DIFF_HIPEAK.out.diff
            if(!skip_peak_annotation){
                BIOC_CHIPPEAKANNO_DIFFHIPEAK(DIFF_HIPEAK.out.diff
                                                        .collect()
                                                        .map{["DiffHiPeak", it]}, gtf)
            }
        }
    }

    emit:
    peak         = ASSIGN_TYPE.out.peak         // channel: [ path(peak) ]
    bedpe        = ASSIGN_TYPE.out.bedpe        // channel: [ path(bedpe) ]
    stats        = stats                        // channel: [ path(stats) ]
    diff         = diff                         // channel: [ path(diff) ]
    versions     = ch_version                   // channel: [ path(version) ]
}
