/*
 * Creating Genomic Features Files for MAPS
 */

include { BIOC_ENZYMECUT               } from '../../../modules/local/bioc/enzyme_cut'
include { MAPS_CUT                     } from '../../../modules/local/maps/cut'
include { MAPS_FEND                    } from '../../../modules/local/maps/fend'
include { SEQLEVELS_STYLE              } from '../../../modules/local/bioc/seqlevelsstyle'
include {
    ENSEMBL_UCSC_CONVERT
        as ENSEMBL_UCSC_CONVERT1;
    ENSEMBL_UCSC_CONVERT
        as ENSEMBL_UCSC_CONVERT2       } from '../../../modules/local/bioc/ensembl_ucsc_convert'
include { UCSC_BIGWIGAVERAGEOVERBED    } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { MAPS_MERGE                   } from '../../../modules/local/maps/merge'
include { MAPS_FEATURE                 } from '../../../modules/local/maps/feature'

workflow MAPS_MULTIENZYME {
    take:
    cool_bin                   // channel: [ val(bin) ]
    genome                     // channel: [ site, [fasta], [chrom_size], [mappability_bw] ]

    main:
    chrom_sizes = genome.map{[ it[2] ]}
    mappability = genome.map{[ it[3] ]}
    if(params.maps_digest_file && params.enzyme.toLowerCase() != "mnase"){
        ch_version = Channel.empty()
        ch_digest = cool_bin.combine(Channel.fromPath(params.maps_digest_file))
        MAPS_FEND(ch_digest, chrom_sizes)
    }else{
        if(params.enzyme.toLowerCase() != "mnase"){
            ch_version = BIOC_ENZYMECUT(
                cool_bin.combine(genome.map{[it[0], it[1]]}), // [bin, site, [fasta]]
                params.enzyme).versions
            ch_digest = BIOC_ENZYMECUT.out.cut
        }else{
            ch_version = MAPS_CUT(
                cool_bin.combine(genome.map{[it[0], it[1]]}), // [bin, site, [fasta]]
                params.enzyme).versions
            ch_digest = MAPS_CUT.out.cut
        }
        MAPS_FEND(ch_digest, chrom_sizes)
    }
    seqlevelsstyle = SEQLEVELS_STYLE(MAPS_FEND.out.bed.map{it[1]}.collect().map{it[0]}).seqlevels_style
    if("$seqlevelsstyle" != "UCSC"){
        ENSEMBL_UCSC_CONVERT1(MAPS_FEND.out.bed)
        ENSEMBL_UCSC_CONVERT2(cool_bin.combine(mappability))
        ch_version = ch_version.mix(ENSEMBL_UCSC_CONVERT1.out.versions)
        UCSC_BIGWIGAVERAGEOVERBED(ENSEMBL_UCSC_CONVERT1.out.tab.map{[['id':'background', 'bin_size':it[0]], it[1]]},
                                    ENSEMBL_UCSC_CONVERT2.out.tab.map{it[1]})
    }else{
        UCSC_BIGWIGAVERAGEOVERBED(MAPS_FEND.out.bed.map{[['id':'background', 'bin_size':it[0]], it[1]]}, mappability)
    }
    ch_version = ch_version.mix(UCSC_BIGWIGAVERAGEOVERBED.out.versions)
    MAPS_MERGE(ch_digest.cross(UCSC_BIGWIGAVERAGEOVERBED.out.tab.map{[it[0].bin_size, it[1]]}).map{[it[0][0], it[0][1], it[1][1]]})

    MAPS_FEATURE(MAPS_MERGE.out.map, chrom_sizes)

    emit:
    bin_feature              = MAPS_FEATURE.out.bin_feature      // channel: [ val(bin_size), path(bin_feature) ]
    versions                 = ch_version                        // channel: [ path(version) ]
}
