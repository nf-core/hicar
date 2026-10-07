/*
 * call compartments by homer
 */

include { HOMER_RUNHICPCA } from '../../../modules/local/homer/run_hic_pca'
include { HOMER_FINDHICCOMPARTMENTS } from '../../../modules/local/homer/find_hic_compartments'
include { BEDGRAPH_TRIM   } from '../../../modules/local/bioc/trimbedgraph'
include {
    UCSC_BEDGRAPHTOBIGWIG
        as UCSC_BEDGRAPHTOBIGWIG_HOMER_COMPARTMENTS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'


workflow HOMER_COMPARTMENTS {
    take:
    tagdir       // channel: [ val(meta), [tagdir] ]
    resolution   // value
    genome       // value

    main:
    HOMER_RUNHICPCA(
        tagdir,
        resolution,
        genome.map{it[0]}.collect()
    )
    ch_version = HOMER_RUNHICPCA.out.versions

    HOMER_FINDHICCOMPARTMENTS(
        HOMER_RUNHICPCA.out.txt,
        genome.map{it[0]}.collect()
    )
    ch_version = ch_version.mix(HOMER_FINDHICCOMPARTMENTS.out.versions)

    BEDGRAPH_TRIM(
        HOMER_FINDHICCOMPARTMENTS.out.bedgraph,
        genome.map{it[1]}.collect()
    )
    ch_version = ch_version.mix(BEDGRAPH_TRIM.out.versions)

    UCSC_BEDGRAPHTOBIGWIG_HOMER_COMPARTMENTS(
        BEDGRAPH_TRIM.out.bedgraph,
        genome.map{it[1]}.collect()
    )
    ch_version = ch_version.mix(UCSC_BEDGRAPHTOBIGWIG_HOMER_COMPARTMENTS.out.versions)

    emit:
    compartments = UCSC_BEDGRAPHTOBIGWIG_HOMER_COMPARTMENTS.out.bigwig  // channel: [ val(meta), path(compartments)]
    versions     = ch_version                                           // channel: [ path(version) ]
}
