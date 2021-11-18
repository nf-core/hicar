/*
 * Createing Stats for mapping results
 */
params.options = [:]

include { SAMTOOLS_SORT                } from '../../modules/nf-core/modules/samtools/sort/main'                     addParams(options: params.options.samtools_sort)
include { SAMTOOLS_INDEX               } from '../../modules/nf-core/modules/samtools/index/main'                    addParams(options: params.options.samtools_index)
include { SAMTOOLS_STATS               } from '../../modules/nf-core/modules/samtools/stats/main'                    addParams(options: params.options.samtools_stats)
include { SAMTOOLS_IDXSTATS            } from '../../modules/nf-core/modules/samtools/idxstats/main'                 addParams(options: params.options.samtools_idxstats)
include { SAMTOOLS_FLAGSTAT            } from '../../modules/nf-core/modules/samtools/flagstat/main'                 addParams(options: params.options.samtools_flagstat)

workflow BAM_STAT {
    take:
    bam          // channel: [ val(meta), path(bam) ]

    main:
    ch_version = SAMTOOLS_SORT(bam).versions
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    SAMTOOLS_STATS(ch_bam_bai, [])
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    SAMTOOLS_IDXSTATS(ch_bam_bai)

    emit:
    stats    = SAMTOOLS_STATS.out.stats           // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat     // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats     // channel: [ val(meta), [ idxstats ] ]
    versions = ch_version                         // channel: [ path(version) ]
}
