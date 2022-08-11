/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHicar.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta,
                            params.gtf, params.bwa_index, params.gene_bed,
                            params.mappability]
for (param in checkPathParamList) {
    if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) {
    ch_input = file("${params.input}", checkIfExists: true)
} else { exit 1, 'Input samplesheet not specified!' }

// set the restriction_sites
def RE_cutsite = [
    "mboi": "^GATC",
    "dpnii": "^GATC",
    "bglii": "^GATCT",
    "hindiii": "^AGCTT",
    "cviqi": "^TAC",
    "msei":"^TAA"]
if (!params.enzyme.toLowerCase() in RE_cutsite){
    exit 1, "Not supported yet!"
}
params.restriction_sites = RE_cutsite[params.enzyme.toLowerCase()]

// if user defined Peaks
ch_anchor_peaks = Channel.empty()
if(params.anchor_peaks){
    ch_anchor_peaks = file("${params.anchor_peaks}", checkIfExists: true)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_circos_config         = file("$projectDir/assets/circos.conf", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_juicer_tools              = file(params.juicer_tools_jar,
                                    checkIfExists: true)
ch_merge_map_py_source       = file(params.merge_map_py_source,
                                    checkIfExists: true)
ch_feature_frag2bin_source   = file(params.feature_frag2bin_source,
                                    checkIfExists: true)
ch_make_maps_runfile_source  = file(params.make_maps_runfile_source,
                                    checkIfExists: true)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { CHECKSUMS } from '../modules/local/checksums'
include { COOLTOOLS_COMPARTMENTS } from '../modules/local/cooltools/eigs-cis'

include { BIOC_CHIPPEAKANNO } from '../modules/local/bioc/chippeakanno'
include { BIOC_CHIPPEAKANNO as BIOC_CHIPPEAKANNO_MAPS } from '../modules/local/bioc/chippeakanno'
include { BIOC_ENRICH } from '../modules/local/bioc/enrich'
include { BIOC_TRACKVIEWER } from '../modules/local/bioc/trackviewer'
include { BIOC_TRACKVIEWER as BIOC_TRACKVIEWER_MAPS } from '../modules/local/bioc/trackviewer'
include { IGV } from '../modules/local/igv'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/preparegenome'
include { BAM_STAT } from '../subworkflows/local/bam_stats'
include { PAIRTOOLS_PAIRE } from '../subworkflows/local/pairtools'
include { COOLER } from '../subworkflows/local/cooler'

include { ATAC_PEAK } from '../subworkflows/local/callatacpeak'
include { TADS } from '../subworkflows/local/tads'
include { COMPARTMENTS } from '../subworkflows/local/compartments'
include { APA } from '../subworkflows/local/aggregate_peak'
include { INTERACTIONS } from '../subworkflows/local/interactions'
include { HIPEAK } from '../subworkflows/local/hipeak'
include { DA } from '../subworkflows/local/differential_analysis'

include { RUN_CIRCOS } from '../subworkflows/local/circos'
include { RUN_CIRCOS as MAPS_CIRCOS } from '../subworkflows/local/circos'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BWA_MEM                     } from '../modules/nf-core/modules/bwa/mem/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { HOMER_MAKETAGDIRECTORY      } from '../modules/nf-core/modules/homer/maketagdirectory/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { SAMTOOLS_MERGE              } from '../modules/nf-core/modules/samtools/merge/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

cool_bin = Channel.fromList(params.cool_bin.tokenize('_'))

workflow HICAR {

    ch_versions         = Channel.empty() // pipeline versions
    ch_circos_files     = Channel.empty() // circos plots
    ch_de_files         = Channel.empty() // Differential analysis
    ch_annotation_files = Channel.empty() // files to be annotated
    ch_multiqc_files    = Channel.empty() // multiQC reports

    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))

    //
    // preprocess: check inputs, check checksums, prepare genome
    //

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_reads = INPUT_CHECK.out.reads

    //
    // check the input fastq files are correct
    //
    CHECKSUMS( ch_reads )
    ch_versions = ch_versions.mix(CHECKSUMS.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: Prepare genome
    //
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.ifEmpty(null))

    //
    // clean reads: run fastqc, trim
    //

    //
    // MODULE: Run FastQC
    //
    if(!params.skip_fastqc){
        FASTQC (
            ch_reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first().ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    }

    //
    // MODULE: trimming
    //
    if(!params.skip_cutadapt){
        CUTADAPT(
            ch_reads
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
        reads4mapping = CUTADAPT.out.reads
    }else{
        reads4mapping = ch_reads
    }

    //
    // prepare paired reads: mapping, pool technique replicates, filter reads, create cooler files
    //

    //
    // MODULE: mapping
    //
    BWA_MEM(
        reads4mapping,
        PREPARE_GENOME.out.bwa_index,
        false
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.ifEmpty(null))

    //
    // Pool the technique replicates
    //
    BWA_MEM.out.bam
                .map{
                    meta, bam ->
                        meta.id = meta.group + "_REP" + meta.replicate
                        [meta.id, meta, bam]
                }
                .groupTuple(by: [0])
                .map{[it[1][0], it[2].flatten()]}
                .set{ mapped_bam }
    //mapped_bam.view()//no branch to multiple and single, need to rename the bam files
    SAMTOOLS_MERGE(mapped_bam, [])
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.ifEmpty(null))

    //
    // MODULE: mapping stats
    //
    BAM_STAT(SAMTOOLS_MERGE.out.bam)
    ch_versions = ch_versions.mix(BAM_STAT.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.idxstats.collect{it[1]}.ifEmpty([]))

    //
    // SUBWORKFLOW: filter reads, output pair (like hic pair), raw (pair), and stats
    //
    PAIRTOOLS_PAIRE(
        SAMTOOLS_MERGE.out.bam,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.digest_genome
    )
    ch_versions = ch_versions.mix(PAIRTOOLS_PAIRE.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(PAIRTOOLS_PAIRE.out.stat.collect().ifEmpty([]))

    //
    // combine bin_size and create cooler file, and dump long_bedpe
    //
    cool_bin.combine(PAIRTOOLS_PAIRE.out.pair)
            .map{bin, meta, pair, px -> [meta, bin, pair, px]}
            .set{cool_input}
    COOLER(
        cool_input,
        PREPARE_GENOME.out.chrom_sizes,
        params.juicer_jvm_params,
        ch_juicer_tools,
        params.long_bedpe_postfix
    )
    ch_versions = ch_versions.mix(COOLER.out.versions.ifEmpty(null))

    //
    // prepare for Homer
    //
    HOMER_MAKETAGDIRECTORY( PAIRTOOLS_PAIRE.out.homerpair, PREPARE_GENOME.out.fasta)

    //
    // 1D peak is required for loops calling
    // calling ATAC peaks, output ATAC narrowPeak and reads in peak
    // or user user predefined peaks
    //
    ATAC_PEAK(// TODO, add TSS heatmap plot, add QC report
        PAIRTOOLS_PAIRE.out.validpair,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gsize,
        PREPARE_GENOME.out.gtf,
        params.method,
        params.anchor_peaks,
        ch_anchor_peaks,
        params.short_bed_postfix
    )
    ch_versions = ch_versions.mix(ATAC_PEAK.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(ATAC_PEAK.out.stats.collect().ifEmpty(null))

    //
    // calling compartments
    //
    if(!params.skip_compartments){
        COMPARTMENTS(
            COOLER.out.mcool,
            params.res_compartments,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(COMPARTMENTS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(COMPARTMENTS.out.mqc.collect().ifEmpty(null))
        ch_circos_files = ch_circos_files.mix(COMPARTMENTS.out.circos)
    }

    //
    // calling TADs
    //
    if(!params.skip_tads){
        TADS(
            COOLER.out.cool,
            params.res_tads,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes
        )
        ch_versions = ch_versions.mix(TADS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(TADS.out.mqc.collect().ifEmpty(null))
        ch_circos_files = ch_circos_files.mix(TADS.out.circos)
    }

    //
    // call interaction loops
    //
    if(!params.skip_interactions){
        INTERACTIONS(
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes,
            PREPARE_GENOME.out.mappability,
            cool_bin,
            PREPARE_GENOME.out.site,
            ATAC_PEAK.out.reads,
            ATAC_PEAK.out.mergedpeak,
            COOLER.out.bedpe,
            HOMER_MAKETAGDIRECTORY.out.tagdir,
            ch_merge_map_py_source,
            ch_feature_frag2bin_source,
            ch_make_maps_runfile_source,
            params.long_bedpe_postfix,
            params.short_bed_postfix,
            params.maps_3d_ext
        )
        ch_versions = ch_versions.mix(INTERACTIONS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(INTERACTIONS.out.mqc.collect().ifEmpty(null))
        ch_annotation_files = ch_annotation_files.mix(INTERACTIONS.out.anno)
        ch_circos_files = ch_circos_files.mix(INTERACTIONS.out.circos)
    }

    //
    // aggregate peak analysis
    //
    if(!params.skip_apa){
        APA(
            COOLER.out.cool,
            COOLER.out.hic,
            ATAC_PEAK.out.mergedpeak
        )
        ch_versions = ch_versions.mix(APA.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(APA.out.mqc.collect().ifEmpty(null))
    }

    //
    // call HiPeak
    // calling high resolution fragments peaks and then call loops
    // this process is time comsuming step
    //
    if(params.high_resolution_R1 && params.method.toLowerCase()=="hicar"){
        HIPEAK(
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.chrom_sizes,
            PREPARE_GENOME.out.digest_genome,
            PREPARE_GENOME.out.mappability,
            PREPARE_GENOME.out.gtf,
            PAIRTOOLS_PAIRE.out.distalpair,
            PAIRTOOLS_PAIRE.out.hdf5,
            ATAC_PEAK.out.peak,
            params.r1_pval_thresh,
            params.short_bed_postfix,
            params.maps_3d_ext
        )
        ch_versions = ch_versions.mix(HIPEAK.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(HIPEAK.out.mqc.ifEmpty(null))
    }

    //
    // Differential analysis
    //
    if(!params.skip_diff_analysis){
        DA(INTERACTIONS.out.loops, COOLER.out.samplebedpe, params.long_bedpe_postfix)
        ch_versions = ch_versions.mix(DA.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(DA.out.mqc.collect().ifEmpty(null))
        ch_annotation_files = ch_annotation_files.mix(DA.out.anno)

      /*  if(ch_diffhicar){
            DIFFHICAR(ch_diffhicar)
            ch_versions = ch_versions.mix(DIFFHICAR.out.versions.ifEmpty(null))
            ch_multiqc_files = ch_multiqc_files.mix(DIFFHICAR.out.stats.collect().ifEmpty([]))
            //annotation
            if(!params.skip_peak_annotation){
                BIOC_CHIPPEAKANNO(DIFFHICAR.out.diff, PREPARE_GENOME.out.gtf)
                ch_versions = ch_versions.mix(BIOC_CHIPPEAKANNO.out.versions.ifEmpty(null))
                if(PREPARE_GENOME.out.ucscname && !params.skip_enrichment){
                    BIOC_ENRICH(
                        BIOC_CHIPPEAKANNO.out.anno.filter{it.size()>0},
                        PREPARE_GENOME.out.ucscname)
                    ch_versions = ch_versions.mix(BIOC_ENRICH.out.versions.ifEmpty(null))
                }
                if(params.virtual_4c){
                    BIOC_CHIPPEAKANNO.out.csv
                        .mix(COOLER.out.mcool
                                .map{meta, mcool -> [meta.bin, mcool]}
                                .groupTuple())
                        .groupTuple()
                        .filter{it.size()>2} //filter the sample without annotation table
                        .map{bin, df -> [bin, df[0], df[1]]}
                        .set{ch_trackviewer}
                    //ch_trackviewer.view()
                    BIOC_TRACKVIEWER(
                        ch_trackviewer,
                        PAIRTOOLS_PAIRE.out.hdf5.collect{it[1]},
                        PREPARE_GENOME.out.gtf,
                        PREPARE_GENOME.out.chrom_sizes,
                        PREPARE_GENOME.out.digest_genome)
                    ch_versions = ch_versions.mix(BIOC_TRACKVIEWER.out.versions.ifEmpty(null))
                }
            }
        }*/
    }

    //
    // Annotation
    //
    if(!params.skip_peak_annotation){
        BIOC_CHIPPEAKANNO(ch_annotation_files, PREPARE_GENOME.out.gtf, params.maps_3d_ext)
        ch_versions = ch_versions.mix(BIOC_CHIPPEAKANNO.out.versions.ifEmpty(null))
    }

    //
    // visualization: virtual_4c
    //

    //
    // visualization: circos
    //
    ch_circos_files.view()
    RUN_CIRCOS(
        ch_circos_files,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.ucscname,
        ch_circos_config
    )
    ch_versions = ch_versions.mix(RUN_CIRCOS.out.versions.ifEmpty(null))

    //
    // visualization: IGV
    //
    INTERACTIONS.out.loops.map{[it[0].id+'.'+it[1]+'.contacts',
                            RelativePublishFolder.getPublishedFolder(workflow,
                                                'MAPS_REFORMAT')+it[2].name]}
        .mix(ATAC_PEAK.out
                    .bws.map{[it[0].id+"_R2",
                        RelativePublishFolder.getPublishedFolder(workflow,
                                            'UCSC_BEDGRAPHTOBIGWIG_PER_GROUP')+it[1].name]})
        .set{ch_trackfiles} // collect track files for igv


/*
    MAPS_CIRCOS(
        INTERACTIONS.out.loops.map{
            meta, bin_size, bedpe ->
                meta.id = "MAPS_PEAK_" + meta.id
                [meta, bedpe]
        },
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.ucscname,
        ch_circos_config
    )
    ch_versions = ch_versions.mix(MAPS_CIRCOS.out.versions.ifEmpty(null))
*/






    //
    // Create igv index.html file
    //
    ch_trackfiles.collect{it.join('\t')}
        .flatten()
        .collectFile(
            name     :'track_files.txt',
            storeDir : params.outdir+'/'+RelativePublishFolder.getPublishedFolder(workflow, 'IGV'),
            newLine  : true, sort:{it[0]})
        .set{ igv_track_files }
    //igv_track_files.view()
    IGV(igv_track_files, PREPARE_GENOME.out.ucscname, RelativePublishFolder.getPublishedFolder(workflow, 'IGV'))

/*
    //
    // Annotate the MAPS peak
    //
    if(!params.skip_peak_annotation){
        INTERACTIONS.out.loops //[]
            .map{meta, bin_size, peak -> [bin_size, peak]}
            .filter{ it[1].readLines().size > 1 }
            .groupTuple()
            .set{ch_maps_anno}
        BIOC_CHIPPEAKANNO_MAPS(ch_maps_anno, PREPARE_GENOME.out.gtf)
        ch_versions = ch_versions.mix(BIOC_CHIPPEAKANNO_MAPS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(BIOC_CHIPPEAKANNO_MAPS.out.png.collect().ifEmpty([]))
        if(params.virtual_4c){
            BIOC_CHIPPEAKANNO_MAPS.out.csv
                .mix(
                    COOLER.out.mcool
                        .map{
                                meta, mcool ->
                                    [meta.bin, mcool]}
                        .groupTuple())
                .groupTuple()
                .filter{it.size()>2} //filter the sample without annotation table
                .map{bin, df -> [bin, df[0], df[1]]}
                .set{ch_maps_trackviewer}
            //ch_maps_trackviewer.view()
            BIOC_TRACKVIEWER_MAPS(
                ch_maps_trackviewer,
                PAIRTOOLS_PAIRE.out.hdf5.collect{it[1]},
                PREPARE_GENOME.out.gtf,
                PREPARE_GENOME.out.chrom_sizes,
                PREPARE_GENOME.out.digest_genome)
            ch_versions = ch_versions.mix(BIOC_TRACKVIEWER_MAPS.out.versions.ifEmpty(null))
        }
    }


*/
    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect().ifEmpty([]))

    if(!params.skip_multiqc){
        //
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowHicar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files
            .flatten()
            .map { it -> if (it) [ it.baseName, it ] }
            .groupTuple()
            .map { it[1][0] }
            .flatten()
            .collect()
            .set { ch_multiqc_files }
        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_versions = ch_versions.mix(MULTIQC.out.versions.ifEmpty(null))
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
