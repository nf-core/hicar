/*
================================================================================
    VALIDATE INPUTS
================================================================================
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
    "cviqi": "^TAC"]
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
================================================================================
    CONFIG FILES
================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml",
                                checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ?
                                Channel.fromPath(params.multiqc_config) :
                                Channel.empty()
ch_circos_config         = file("$projectDir/assets/circos.conf",
                                checkIfExists: true)

/*
================================================================================
    TOOLS SOURCE FILE
================================================================================
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
================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Local to the pipeline
//
include { CHECKSUMS } from '../modules/local/checksums'
include { DIFFHICAR } from '../modules/local/bioc/diffhicar'
include { BIOC_CHIPPEAKANNO } from '../modules/local/bioc/chippeakanno'
include { BIOC_CHIPPEAKANNO as BIOC_CHIPPEAKANNO_MAPS } from '../modules/local/bioc/chippeakanno'
include { BIOC_ENRICH } from '../modules/local/bioc/enrich'
include { BIOC_TRACKVIEWER } from '../modules/local/bioc/trackviewer'
include { BIOC_TRACKVIEWER as BIOC_TRACKVIEWER_MAPS } from '../modules/local/bioc/trackviewer'
include { IGV } from '../modules/local/igv'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME } from '../subworkflows/local/preparegenome'
include { BAM_STAT } from '../subworkflows/local/bam_stats'
include { PAIRTOOLS_PAIRE } from '../subworkflows/local/pairtools'
include { COOLER } from '../subworkflows/local/cooler'
include { ATAC_PEAK } from '../subworkflows/local/callatacpeak'
include { R1_PEAK } from '../subworkflows/local/calldistalpeak'
include { HI_PEAK } from '../subworkflows/local/hipeak'
include { MAPS_MULTIENZYME } from '../subworkflows/local/multienzyme'
include { MAPS_PEAK } from '../subworkflows/local/maps_peak'
include { RUN_CIRCOS } from '../subworkflows/local/circos'
include { RUN_CIRCOS as MAPS_CIRCOS } from '../subworkflows/local/circos'
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { BWA_MEM                     } from '../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_MERGE              } from '../modules/nf-core/modules/samtools/merge/main'

/*
================================================================================
    RUN MAIN WORKFLOW
================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

cool_bin = Channel.fromList(params.cool_bin.tokenize('_'))

workflow HICAR {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_reads = INPUT_CHECK.out.reads

    //
    // check the input fastq files are correct and produce checksum for GEO submission
    //
    CHECKSUMS( ch_reads )
    ch_versions = ch_versions.mix(CHECKSUMS.out.versions.ifEmpty(null))

    //
    // SUBWORKFLOW: Prepare genome
    //
    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.ifEmpty(null))

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
        ch_juicer_tools
    )
    ch_versions = ch_versions.mix(COOLER.out.versions.ifEmpty(null))

    //
    // calling ATAC peaks, output ATAC narrowPeak and reads in peak
    // or user user predefined peaks
    //
    ATAC_PEAK(
        PAIRTOOLS_PAIRE.out.validpair,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gsize,
        PREPARE_GENOME.out.gtf,
        params.method,
        params.anchor_peaks,
        ch_anchor_peaks
    )
    ch_versions = ch_versions.mix(ATAC_PEAK.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(ATAC_PEAK.out.stats.collect().ifEmpty([]))

    //
    // calling distal peaks: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]
    //
    background = MAPS_MULTIENZYME(PREPARE_GENOME.out.fasta,
                                    cool_bin,
                                    PREPARE_GENOME.out.chrom_sizes,
                                    ch_merge_map_py_source,
                                    ch_feature_frag2bin_source).bin_feature
    ch_versions = ch_versions.mix(MAPS_MULTIENZYME.out.versions.ifEmpty(null))
    reads_peak   = ATAC_PEAK.out.reads
                            .map{ meta, reads ->
                                    [meta.id, reads]} // here id is group
                            .combine(ATAC_PEAK.out.mergedpeak)// group, reads, peaks
                            .cross(COOLER.out.bedpe.map{[it[0].id, it[0].bin, it[1]]})// group, bin, bedpe
                            .map{ short_bed, long_bedpe -> //[bin_size, group, macs2, long_bedpe, short_bed]
                                    [long_bedpe[1], short_bed[0], short_bed[2], long_bedpe[2], short_bed[1]]}
    background.cross(reads_peak)
                .map{ background, reads -> //[group, bin_size, macs2, long_bedpe, short_bed, background]
                        [[id:reads[1]], background[0], reads[2], reads[3], reads[4], background[1]]}
                .set{ maps_input }
    MAPS_PEAK(
        maps_input,
        ch_make_maps_runfile_source,
        PREPARE_GENOME.out.chrom_sizes,
        params.juicer_jvm_params,
        ch_juicer_tools)
    ch_versions = ch_versions.mix(MAPS_PEAK.out.versions.ifEmpty(null))
    ch_multiqc_files = ch_multiqc_files.mix(MAPS_PEAK.out.stats.collect().ifEmpty([]))

    MAPS_CIRCOS(
        MAPS_PEAK.out.peak.map{
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

    MAPS_PEAK.out.peak.map{[it[0].id+'.'+it[1]+'.contacts',
                            RelativePublishFolder.getPublishedFolder(workflow,
                                                'MAPS_REFORMAT')+it[2].name]}
        .mix(ATAC_PEAK.out
                    .bws.map{[it[0].id+"_R2",
                        RelativePublishFolder.getPublishedFolder(workflow,
                                            'UCSC_BEDGRAPHTOBIGWIG_PER_GROUP')+it[1].name]})
        .set{ch_trackfiles} // collect track files for igv

    //
    // calling R1 peaks, output R1 narrowPeak and reads in peak
    //
    if(params.high_resolution_R1 && params.method.toLowerCase()=="hicar"){
        R1_PEAK(
            PAIRTOOLS_PAIRE.out.distalpair,
            PREPARE_GENOME.out.chrom_sizes,
            PREPARE_GENOME.out.digest_genome,
            PREPARE_GENOME.out.gtf,
            params.r1_pval_thresh
        )
        ch_versions = ch_versions.mix(R1_PEAK.out.versions.ifEmpty(null))
        ch_trackfiles = ch_trackfiles.mix(
            R1_PEAK.out.bws.map{[it[0].id+"_R1",
                RelativePublishFolder.getPublishedFolder(workflow,
                    'UCSC_BEDGRAPHTOBIGWIG_PER_R1_GROUP')+it[1].name]})

        // merge ATAC_PEAK with R1_PEAK by group id
        distalpair = PAIRTOOLS_PAIRE.out.hdf5.map{meta, bed -> [meta.group, bed]}
                                            .groupTuple()
        grouped_reads_peak = ATAC_PEAK.out.peak.map{[it[0].id, it[1]]}
                                .join(R1_PEAK.out.peak.map{[it[0].id, it[1]]})
                                .join(distalpair)
                                .map{[[id:it[0]], it[1], it[2], it[3]]}
        HI_PEAK(
            grouped_reads_peak,
            PREPARE_GENOME.out.chrom_sizes,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.digest_genome,
            MAPS_MULTIENZYME.out.mappability,
            params.skip_peak_annotation,
            params.skip_diff_analysis
        )
        ch_versions = ch_versions.mix(HI_PEAK.out.versions.ifEmpty(null))
        ch_trackfiles = ch_trackfiles.mix(
            HI_PEAK.out.bedpe
                    .map{[it[0].id+"_HiPeak",
                        RelativePublishFolder.getPublishedFolder(workflow,
                                            'ASSIGN_TYPE')+it[1].name]})

        RUN_CIRCOS(
            HI_PEAK.out.bedpe,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.chrom_sizes,
            PREPARE_GENOME.out.ucscname,
            ch_circos_config
        )
        ch_versions = ch_versions.mix(RUN_CIRCOS.out.versions.ifEmpty(null))
    }

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

    //
    // Annotate the MAPS peak
    //
    if(!params.skip_peak_annotation){
        MAPS_PEAK.out.peak //[]
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
                .filter{it.size()>2} //filter to remove the groups with no different results.
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

    //
    // Differential analysis
    //
    if(!params.skip_diff_analysis){
        MAPS_PEAK.out.peak //[]
            .map{meta, bin_size, peak -> [bin_size, peak]}
            .groupTuple()
            .cross(COOLER.out.samplebedpe.map{[it[0].bin, it[1]]}.groupTuple())
            .map{ peak, long_bedpe ->
                [peak[0], peak[1].flatten(), long_bedpe[1].flatten()] }//bin_size, meta, peak, long_bedpe
            .groupTuple()
            .map{[it[0], it[1].flatten().unique(), it[2].flatten()]}
            .filter{it[1].size > 1} // filter by the bedpe files. Single bedpe means single group, no need to do differential analysis
            .set{ch_diffhicar}
        //ch_diffhicar.view()
        if(ch_diffhicar){
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
                        .filter{it.size()>2} //filter to remove the groups with no different results.
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
        }
    }

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
================================================================================
    COMPLETION EMAIL AND SUMMARY
================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
================================================================================
    THE END
================================================================================
*/
