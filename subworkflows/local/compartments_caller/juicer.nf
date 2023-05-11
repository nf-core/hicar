/*
 * call compartments by Juicer
 */

include { JUICER_EIGENVECTOR      } from '../../../modules/local/juicer/eigenvector'
include { BEDGRAPH_TRIM
    as WIG_TRIM } from '../../../modules/local/bioc/trimbedgraph'
include { UCSC_BEDGRAPHTOBIGWIG
    as UCSC_BEDGRAPHTOBIGWIG_JUICER_EIGENVECTOR   } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'


workflow JUICER_COMPARTMENTS {
    take:
    matrix                 // channel: [ val(meta), [cool] ]
    resolution             // value
    additional_param       // value: [ [juicer_tools_jar], [chrom_size] ]

    main:
    JUICER_EIGENVECTOR(
        matrix.combine(additional_param),
        resolution
    )
    ch_version = JUICER_EIGENVECTOR.out.versions

    WIG_TRIM(
        JUICER_EIGENVECTOR.out.wig,
        additional_param.map{it[1]}.collect()
    )
    ch_version = ch_version.mix(WIG_TRIM.out.versions)

    UCSC_BEDGRAPHTOBIGWIG_JUICER_EIGENVECTOR(
        WIG_TRIM.out.bedgraph,
        additional_param.map{it[1]}.collect()
    )
    ch_compartments = UCSC_BEDGRAPHTOBIGWIG_JUICER_EIGENVECTOR.out.bigwig
        .map{
            meta,bw ->
                smeta = meta.clone()
                smeta.id = smeta.id.replaceAll("_${meta.bin}_${resolution}", "")
                [smeta, bw]

        }
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_JUICER_EIGENVECTOR.out.versions)

    emit:
    eigenvectors   = JUICER_EIGENVECTOR.out.compartments         // channel: [ val(meta), path(files)]
    compartments   = ch_compartments                             // channel: [ val(meta), path(bigwig)]
    versions  = ch_version                                       // channel: [ path(version) ]
}
