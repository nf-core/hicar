/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowHicar.initialise(params, log)

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
    PIPELINE CONTROLER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// check the tool is used
def checkToolsUsedInDownstream(tool, params){
    return (
        (params.interactions_tool == tool && !params.skip_interactions) ||
        (params.tad_tool == tool && !params.skip_tads) ||
        (params.compartments_tool == tool && !params.skip_compartments) ||
        (params.apa_tool == tool && params.do_apa) ||
        (params.da_tool == tool && !params.skip_diff_analysis) ||
        (params.v4c_tool == tool && params.create_virtual_4c) ||
        (params.tfea_tool == tool && params.do_tfea)
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
// custom config file
ch_circos_config         = file("$projectDir/assets/circos.conf", checkIfExists: true)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_juicer_tools              = Channel.fromPath(params.juicer_tools_jar,
                                    checkIfExists: true).collect()
ch_hic_tools                 = Channel.fromPath(params.hic_tools_jar,
                                    checkIfExists: true).collect()
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local to the pipeline
//
include { CHECKSUMS                                } from '../modules/local/checksums'
include { HOMER_INSTALL                            } from '../modules/local/homer/install'
include { RECENTER_PEAK                            } from '../modules/local/hicexplorer/recenterpeak'
include { HICEXPLORER_CHICQUALITYCONTROL           } from '../modules/local/hicexplorer/chicqualitycontrol'
include { HICEXPLORER_CHICVIEWPOINTBACKGROUNDMODEL } from '../modules/local/hicexplorer/chicviewpointbackgroundmodel'
include { HICEXPLORER_CHICVIEWPOINT                } from '../modules/local/hicexplorer/chicviewpoint'
include { JUICER_ADDNORM                           } from '../modules/local/juicer/addnorm'
include { BIOC_CHIPPEAKANNO                        } from '../modules/local/bioc/chippeakanno'
include { BIOC_ENRICH                              } from '../modules/local/bioc/enrich'
include { IGV                                      } from '../modules/local/igv'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME  } from '../subworkflows/local/preparegenome'
include { BAM_STAT        } from '../subworkflows/local/bam_stats'
include { PAIRTOOLS_PAIRE } from '../subworkflows/local/pairtools'
include { COOLER          } from '../subworkflows/local/cooler'

include { ATAC_PEAK       } from '../subworkflows/local/callatacpeak'
include { TADS            } from '../subworkflows/local/tads'
include { COMPARTMENTS    } from '../subworkflows/local/compartments'
include { APA             } from '../subworkflows/local/apa'
include { INTERACTIONS    } from '../subworkflows/local/interactions'
include { HIPEAK          } from '../subworkflows/local/hipeak'
include { DA              } from '../subworkflows/local/differential_analysis'
include { V4C             } from '../subworkflows/local/v4c'
include { TFEA            } from '../subworkflows/local/tfea'

include { RUN_CIRCOS      } from '../subworkflows/local/circos'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BWA_MEM                     } from '../modules/nf-core/bwa/mem/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CUTADAPT                    } from '../modules/nf-core/cutadapt/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { CAT_CAT                     } from '../modules/nf-core/cat/cat/main'
include { HOMER_MAKETAGDIRECTORY      } from '../modules/nf-core/homer/maketagdirectory/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_MERGE              } from '../modules/nf-core/samtools/merge/main'

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
    ch_multiqc_files    = Channel.empty() // multiQC reports
    ch_circos_files     = Channel.empty() // circos plots
    ch_de_files         = Channel.empty() // Differential analysis
    ch_annotation_files = Channel.empty() // files to be annotated
    ch_v4c_files        = Channel.empty() // files used for v4c events
    ch_loop_1d_peak     = Channel.empty() // 1D peak for interactions/loops calling
    ch_loop_matrix      = Channel.empty() // matrix for interactions/loops calling
    ch_apa_matrix       = Channel.empty() // matrix for apa analysis
    ch_tad_matrix       = Channel.empty() // matrix for tad analysis
    ch_comp_matrix      = Channel.empty() // matrix for compartment calling
    ch_tfea_bed         = Channel.empty() // peaks for TFEA
    ch_loop_additional  = Channel.empty() // additional inputs for interaction/loops calling
    ch_apa_additional   = Channel.empty() // additional inputs for apa
    ch_tad_additional   = Channel.empty() // additional inputs for tad
    ch_comp_additional  = Channel.empty() // additional inputs for A/B compartments
    ch_tfea_additional  = Channel.empty() // additional inputs for TFEA

    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_config)
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))

    //
    // preprocess: check inputs, check checksums, prepare genome
    //

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_reads = Channel.fromSamplesheet("input")
        .map{
            meta, fastq_1, fastq_2 ->
                meta_copy = meta - meta.subMap(['id', 'single_end', 'group']) + [id: meta.group.toString().replaceAll("\\.", "_") + "_REP" + meta.replicate + "_T" + meta.techniquereplicate, single_end: false, group: meta.group.toString().replaceAll("\\.", "_")]
                [meta_copy, [fastq_1, fastq_2]]
        }

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
    SAMTOOLS_MERGE(mapped_bam, [], [])
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
        PREPARE_GENOME.out.digest_genome,
        params.resample_pairs
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
        ch_hic_tools,
        params.long_bedpe_postfix
    )
    ch_versions = ch_versions.mix(COOLER.out.versions.ifEmpty(null))

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
    ch_multiqc_files = ch_multiqc_files.mix(ATAC_PEAK.out.stats.collect().ifEmpty([]))
    ch_tfea_bed = ATAC_PEAK.out.peak.map{[it[0], 'R2', it[1]]}

    //
    // prepare for MAPS
    //
    if(params.interactions_tool == 'maps'){
        ch_loop_matrix = COOLER.out.bedpe
        ch_loop_1d_peak = ATAC_PEAK.out.reads.combine(ATAC_PEAK.out.mergedpeak)
        ch_loop_additional = PREPARE_GENOME.out.site
                                .combine(PREPARE_GENOME.out.fasta)
                                .combine(PREPARE_GENOME.out.chrom_sizes)
                                .combine(PREPARE_GENOME.out.mappability)
                                .collect()
    }

    //
    // prepare for HiC-DC+
    //
    if(params.interactions_tool == 'hicdcplus'){
        ch_loop_matrix = COOLER.out.bedpe
        ch_loop_1d_peak = ATAC_PEAK.out.reads
        ch_loop_additional = PREPARE_GENOME.out.site
                                .combine(PREPARE_GENOME.out.fasta)
                                .combine(PREPARE_GENOME.out.chrom_sizes)
                                .combine(PREPARE_GENOME.out.mappability)
    }

    //
    // prepare for peakachu
    //
    if(params.interactions_tool == 'peakachu' && params.method == 'HiCAR'){
        ch_loop_matrix = COOLER.out.cool
    }

    //
    // prepare for HiCExplorer
    //
    if(checkToolsUsedInDownstream('hicexplorer', params)){
        if(params.compartments_tool == 'hicexplorer'){
            ch_comp_matrix = COOLER.out.cool
            ch_comp_additional = PREPARE_GENOME.out.chrom_sizes
        }
        if(params.tad_tool == 'hicexplorer'){
            ch_tad_matrix = COOLER.out.cool
            ch_tad_additional = PREPARE_GENOME.out.chrom_sizes
        }
        if(params.apa_tool == 'hicexplorer'){
            ch_apa_matrix = COOLER.out.cool
        }
        if((params.da_tool == 'hicexplorer' && !params.skip_diff_analysis)||
            (params.v4c_tool == 'hicexplorer' && params.create_virtual_4c)){
            // get viewpoint, input is the merged peaks, [bed]
            RECENTER_PEAK(ATAC_PEAK.out.mergedpeak)

            // create background model, input is [bin_size, [cool], [recentered_peak]]
            COOLER.out.cool
                .map{ meta, cool ->
                            [meta.bin, cool]}
                .groupTuple()
                .combine(RECENTER_PEAK.out.peak)
                .set{ch_chic_explorer}
            HICEXPLORER_CHICQUALITYCONTROL(ch_chic_explorer)
            HICEXPLORER_CHICVIEWPOINTBACKGROUNDMODEL(HICEXPLORER_CHICQUALITYCONTROL.out.referencepoints)

            // create interaction file, input is [bin_size, [cool], [recentered_peak], [background_model]]
            HICEXPLORER_CHICVIEWPOINT(HICEXPLORER_CHICQUALITYCONTROL.out.referencepoints.combine(HICEXPLORER_CHICVIEWPOINTBACKGROUNDMODEL.out.background_model, by: 0))
        }
    }

    //
    // prepare for cooltools
    //
    if(checkToolsUsedInDownstream('cooltools', params)){
        if(params.compartments_tool == 'cooltools'){
            ch_comp_matrix = COOLER.out.mcool // cooltools ask the resolution match the tiled genome
            ch_comp_additional = PREPARE_GENOME.out.fasta.combine(PREPARE_GENOME.out.chrom_sizes)
        }
        if(params.tad_tool == 'cooltools'){
            ch_tad_matrix = COOLER.out.cool
        }
        if(params.apa_tool == 'cooltools'){
            ch_apa_matrix = COOLER.out.cool
        }
    }

    //
    // prepare for Homer
    //
    if(checkToolsUsedInDownstream('homer', params)){
        if(!workflow.containerEngine){//homer data folder writable
            HOMER_INSTALL(
                PREPARE_GENOME.out.ucscname
            )
            homer_genome = PREPARE_GENOME.out.ucscname // values: eg. 'hg38'
            homer_done = HOMER_INSTALL.out.output
        }else{// using custom genomes and annotation files 'on-the-fly'
            homer_genome = PREPARE_GENOME.out.fasta
            homer_done = true
        }

        if(homer_done){// force wait Homer install done
            CAT_CAT(
                PAIRTOOLS_PAIRE.out.homerpair.map{
                    meta, pairs ->
                    [[id:meta.group], pairs]
                }.groupTuple(by:0)
            )
            HOMER_MAKETAGDIRECTORY(
                CAT_CAT.out.file_out,
                PREPARE_GENOME.out.fasta
            )
            if(params.tad_tool == 'homer'){
                ch_tad_matrix = HOMER_MAKETAGDIRECTORY.out.tagdir
                ch_tad_additional = homer_genome
            }
            if(params.compartments_tool == 'homer'){
                ch_comp_matrix = HOMER_MAKETAGDIRECTORY.out.tagdir
                ch_comp_additional = homer_genome.combine(PREPARE_GENOME.out.chrom_sizes)
            }
            if(params.tfea_tool == 'homer'){
                ch_tfea_additional = homer_genome
            }
            ch_versions = ch_versions.mix(CAT_CAT.out.versions.ifEmpty(null))
            ch_versions = ch_versions.mix(HOMER_MAKETAGDIRECTORY.out.versions.ifEmpty(null))
        }
    }

    //
    // prepare for JuicerBox
    //
    if(checkToolsUsedInDownstream('juicebox', params)){
        ch_norm_hic = JUICER_ADDNORM(COOLER.out.hic, ch_hic_tools).hic
        if(params.compartments_tool == 'juicebox'){
            ch_comp_matrix = ch_norm_hic
            ch_comp_additional  = ch_hic_tools.combine(PREPARE_GENOME.out.chrom_sizes)
        }
        if(params.apa_tool == 'juicebox'){
            ch_apa_matrix = ch_norm_hic
            ch_apa_additional = ch_juicer_tools
        }
    }

    //
    // calling compartments
    //
    if(!params.skip_compartments){
        COMPARTMENTS(
            ch_comp_matrix,
            params.res_compartments,
            ch_comp_additional,
            cool_bin
        )
        ch_versions = ch_versions.mix(COMPARTMENTS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(COMPARTMENTS.out.mqc.collect().ifEmpty([]))
        ch_circos_files = ch_circos_files.mix(COMPARTMENTS.out.circos)
    }

    //
    // calling TADs
    //
    if(!params.skip_tads){
        TADS(
            ch_tad_matrix,
            params.res_tads,
            ch_tad_additional,
            cool_bin
        )
        ch_versions = ch_versions.mix(TADS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(TADS.out.mqc.collect().ifEmpty([]))
        ch_circos_files = ch_circos_files.mix(TADS.out.circos)
    }

    //
    // call interaction loops
    //
    if(!params.skip_interactions){
        INTERACTIONS(
            ch_loop_matrix,
            ch_loop_1d_peak,
            ch_loop_additional
        )
        ch_versions = ch_versions.mix(INTERACTIONS.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(INTERACTIONS.out.mqc.collect().ifEmpty([]))
        ch_annotation_files = ch_annotation_files.mix(INTERACTIONS.out.anno.map{[params.interactions_tool+'/'+it[0], it[1]]})
        ch_circos_files = ch_circos_files.mix(INTERACTIONS.out.circos)
        ch_de_files = ch_de_files.mix(INTERACTIONS.out.loops)
    }

    //
    // aggregate peak analysis
    //
    if(params.do_apa){
        ch_apa_peak = params.apa_peak ? Channel.fromPath( params.apa_peak, checkIfExists: true ) : ATAC_PEAK.out.mergedpeak
        if(params.apa_tool=='juicebox'){
            ch_apa_additional = INTERACTIONS.out.mergedloops.combine(ch_apa_additional).unique()
        }
        APA(
            ch_apa_matrix,
            ch_apa_peak,
            ch_apa_additional
        )
        ch_versions = ch_versions.mix(APA.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(APA.out.mqc.collect().ifEmpty([]))
    }

    //
    // call HiPeak
    // calling high resolution fragments peaks and then call loops
    // this process is time comsuming step
    //
    if(params.call_high_peak && params.method.toLowerCase()=="hicar"){
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
        ch_multiqc_files = ch_multiqc_files.mix(HIPEAK.out.mqc.collect().ifEmpty([]))
        ch_tfea_bed = ch_tfea_bed.mix(HIPEAK.out.fragmentPeak.map{[it[0], 'R1', it[1]]})
        ch_tfea_bed = ch_tfea_bed.mix(HIPEAK.out.bed4tfea)
    }

    //
    // Differential analysis
    //
    if(!params.skip_diff_analysis){
        if(params.da_tool=='hicexplorer'){
            DA(
                HICEXPLORER_CHICVIEWPOINT.out.interactions,
                []
            )
        }else{
            DA(
                ch_de_files,
                COOLER.out.samplebedpe
            )
        }
        ch_versions = ch_versions.mix(DA.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(DA.out.mqc.collect().ifEmpty([]))
        ch_annotation_files = ch_annotation_files.mix(DA.out.anno.map{[params.da_tool+'/'+it[1], it[2]]})
        ch_v4c_files = ch_v4c_files.mix(DA.out.anno.map{[it[0], it[2]]})
    }

    //
    // Motif analysis: for R2 reads and R1 reads
    //
    if(params.do_tfea){
        if(params.tfea_tool=='atacseqtfea'){
            ch_tfea_bed = ch_tfea_bed.map{[[id:it[1]], it[2]]}.groupTuple(by:0)
                            .combine(
                                SAMTOOLS_MERGE.out.bam
                                        .map{[it[1]]}
                                        .collect().toList())
            ch_tfea_additional = PREPARE_GENOME.out.ucscname
                                    .combine(PREPARE_GENOME.out.fasta)
                                    .combine(PREPARE_GENOME.out.gtf)
        }
        TFEA(ch_tfea_bed, ch_tfea_additional)
        ch_versions = ch_versions.mix(TFEA.out.versions.ifEmpty(null))
    }

    //
    // Annotation
    //
    if(!params.skip_peak_annotation){
        BIOC_CHIPPEAKANNO(ch_annotation_files.combine(ATAC_PEAK.out.mergedpeak), PREPARE_GENOME.out.gtf, params.maps_3d_ext)
        ch_versions = ch_versions.mix(BIOC_CHIPPEAKANNO.out.versions.ifEmpty(null))
    }

    //
    // visualization: virtual_4c
    //
    if(params.create_virtual_4c){
        if(params.v4c_tool == 'hicexplorer'){
            if(params.da_tool == 'hicexplorer'){
                V4C(
                    DA.out.diff,
                    [],[],[],[]
                )
                ch_versions = ch_versions.mix(V4C.out.versions.ifEmpty(null))
            }
        }else{
            if(params.v4c_tool=='trackviewer'){
                ch_cooler_for_v4c = COOLER.out.mcool
            }else{
                ch_cooler_for_v4c = COOLER.out.cool
            }
            ch_v4c_files
                .mix(
                    ch_cooler_for_v4c
                        .map{
                                meta, cool ->
                                    [meta.bin, cool]}
                        .groupTuple())
                .groupTuple()
                .filter{it[1].size()>1} // in case the ch_v4c files is empty
                .map{bin, df -> [bin, df[0], df[1]]} // [bin, cool, bedpe]
                .combine(ATAC_PEAK.out.mergedpeak)
                .set{ch_v4c}
            V4C(
                ch_v4c,
                PAIRTOOLS_PAIRE.out.hdf5.collect{it[1]},
                PREPARE_GENOME.out.gtf,
                PREPARE_GENOME.out.chrom_sizes,
                PREPARE_GENOME.out.digest_genome)
            ch_versions = ch_versions.mix(V4C.out.versions.ifEmpty(null))
        }
    }

    //
    // visualization: circos
    //
    RUN_CIRCOS(
        ch_circos_files.groupTuple(),
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.ucscname,
        ch_circos_config
    )
    ch_versions = ch_versions.mix(RUN_CIRCOS.out.versions.ifEmpty(null))

    //
    // visualization: IGV, Create igv index.html file
    //
    def bedpe_module_name = 'MAPS_REFORMAT'
    switch(params.interactions_tool){
        case 'maps':
            bedpe_module_name = 'MAPS_REFORMAT'
            break
        case 'hicdcplus':
            bedpe_module_name = 'HICDCPLUS_CALL_LOOPS'
            break
        case 'peakachu':
            bedpe_module_name = 'PEAKACHU_SCORE'
            break
    }
    INTERACTIONS.out.loops.map{[it[0].id+'.'+it[1]+'.contacts',
                            RelativePublishFolder.getPublishedFolder(workflow,
                                                bedpe_module_name)+it[2].name]}
        .mix(ATAC_PEAK.out
                    .bws.map{[it[0].id+"_R2",
                        RelativePublishFolder.getPublishedFolder(workflow,
                                            'UCSC_BEDGRAPHTOBIGWIG_PER_GROUP')+it[1].name]})
        .set{ch_trackfiles} // collect track files for igv

    ch_trackfiles.collect{it.join('\t')}
        .flatten()
        .collectFile(
            name     :'track_files.txt',
            storeDir : params.outdir+'/'+RelativePublishFolder.getPublishedFolder(workflow, 'IGV'),
            newLine  : true, sort:{it[0]})
        .set{ igv_track_files }
    IGV(igv_track_files, PREPARE_GENOME.out.ucscname, RelativePublishFolder.getPublishedFolder(workflow, 'IGV'))

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.yml.collect().ifEmpty([]))

    ch_multiqc_files
        .flatten()
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_multiqc_files }

    if(!params.skip_multiqc){
        //
        // MODULE: MultiQC
        //

        workflow_summary    = WorkflowHicar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowHicar.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
