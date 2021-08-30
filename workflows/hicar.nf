/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowHicar.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath("${params.input}").splitCsv(header: true, sep:",") } else { exit 1, 'Input samplesheet not specified!' }

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

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// Extract parameters from params.modules
def getParam(modules, module) {
    return modules[module]?:[:]
}
def getSubWorkFlowParam(modules, mods) {
    def Map options = [:]
    mods.each{
        val ->
        options[val] = modules[val]?:[:]
    }
    return options
}

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS  } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { CHECKSUMS              } from '../modules/local/checksums' addParams( options: getParam(modules, 'checksums') )
include { DIFFHICAR              } from '../modules/local/bioc/diffhicar' addParams(options: getParam(modules, 'diffhicar'))
include { BIOC_CHIPPEAKANNO      } from '../modules/local/bioc/chippeakanno' addParams(options: getParam(modules, 'chippeakanno'))
include { BIOC_CHIPPEAKANNO
    as BIOC_CHIPPEAKANNO_MAPS    } from '../modules/local/bioc/chippeakanno' addParams(options: getParam(modules, 'chippeakanno_maps'))
include { BIOC_ENRICH            } from '../modules/local/bioc/enrich' addParams(options: getParam(modules, 'enrichment'))
include { BIOC_TRACKVIEWER       } from '../modules/local/bioc/trackviewer' addParams(options: getParam(modules, 'trackviewer'))
include { BIOC_TRACKVIEWER
    as BIOC_TRACKVIEWER_MAPS     } from '../modules/local/bioc/trackviewer' addParams(options: getParam(modules, 'trackviewer_maps'))

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME         } from '../subworkflows/local/preparegenome' addParams ( options: getSubWorkFlowParam(modules, ['gunzip', 'gtf2bed', 'chromsizes', 'genomefilter', 'bwa_index', 'gffread', 'digest_genome']) )
include { BAM_STAT               } from '../subworkflows/local/bam_stats' addParams(options: getSubWorkFlowParam(modules, ['samtools_sort', 'samtools_index', 'samtools_stats', 'samtools_flagstat', 'samtools_idxstats']))
include { PAIRTOOLS_PAIRE        } from '../subworkflows/local/pairtools' addParams(options: getSubWorkFlowParam(modules, ['paritools_dedup', 'pairtools_flip', 'pairtools_parse', 'pairtools_restrict', 'pairtools_select', 'pairtools_select_long', 'pairtools_sort', 'pairix', 'reads_stat', 'reads_summary', 'pairsqc', 'pairsplot']))
include { COOLER                 } from '../subworkflows/local/cooler' addParams(options: getSubWorkFlowParam(modules, ['cooler_cload', 'cooler_merge', 'cooler_zoomify', 'cooler_dump_per_group', 'cooler_dump_per_sample', 'dumpintrareads_per_group', 'dumpintrareads_per_sample']))
include { ATAC_PEAK              } from '../subworkflows/local/callatacpeak' addParams(options: getSubWorkFlowParam(modules, ['pairtools_select', 'pairtools_select_short', 'merge_reads', 'shift_reads', 'macs2_atac', 'dump_reads_per_group', 'dump_reads_per_sample', 'merge_peak', 'bedtools_genomecov_per_group', 'bedtools_genomecov_per_sample', 'ucsc_bedgraphtobigwig_per_group', 'ucsc_bedgraphtobigwig_per_sample']))
include { MAPS_MULTIENZYME       } from '../subworkflows/local/multienzyme'   addParams(options: getSubWorkFlowParam(modules, ['maps_cut', 'maps_fend', 'genmap_index', 'genmap_mappability', 'ucsc_wigtobigwig', 'maps_mapability', 'maps_merge', 'maps_feature', 'ensembl_ucsc_convert']))
include { MAPS_PEAK              } from '../subworkflows/local/maps_peak' addParams(options: getSubWorkFlowParam(modules, ['maps_maps', 'maps_callpeak', 'maps_reformat']))

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC         } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { CUTADAPT       } from '../modules/nf-core/modules/cutadapt/main' addParams(options: getParam(modules, 'cutadapt'))
include { BWA_MEM        } from '../modules/nf-core/modules/bwa/mem/main'  addParams(options: getParam(modules, 'bwa_mem'))
include { SAMTOOLS_MERGE } from '../modules/nf-core/modules/samtools/merge/main'  addParams(options: getParam(modules, 'samtools_merge'))
include { MULTIQC        } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

// Parse input
ch_fastq = ch_input.map{
    row ->
        if(!row.group) { exit 1, 'Input samplesheet must contain 'group' column!' }
        if(!row.replicate) { exit 1, 'Input samplesheet must contain 'replicate' column!' }
        if(!row.fastq_1) { exit 1, 'Input samplesheet must contain 'fastq_1' column!' }
        if(!row.fastq_2) { exit 1, 'Input samplesheet must contain 'fastq_2' column!' }
        if(row.id) { exit 1, 'Input samplesheet can not contain 'id' column!' }
        fastq1 = file(row.remove("fastq_1"), checkIfExists: true)
        fastq2 = file(row.remove("fastq_2"), checkIfExists: true)
        meta = row
        meta.id = row.group + "_REP" + row.replicate
        [meta.id, meta, [fastq1, fastq2]]
}
// rename the input if there are technique duplicates
ch_fastq.groupTuple(by:[0])
        .map{
            id, meta, fq ->
                meta.eachWithIndex{
                    entry, index ->
                        entry.id = entry.id + "_T" + index
                        entry
            }
            [id, meta, fq]
        }.transpose()
        .map{[it[1], it[2]]}
        .set{ ch_reads }

//ch_reads.view()
cool_bin = Channel.fromList(params.cool_bin.tokenize('_'))

workflow HICAR {

    ch_software_versions = Channel.empty()

    //
    // check the input fastq files are correct and produce checksum for GEO submission
    //
    CHECKSUMS( ch_reads )

    //
    // SUBWORKFLOW: Prepare genome
    //
    PREPARE_GENOME()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.version.ifEmpty(null))

    //
    // MODULE: Run FastQC
    //
    if(!params.skip_fastqc){
        FASTQC (
            ch_reads
        )
        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    }

    //
    // MODULE: trimming
    //
    CUTADAPT(
        ch_reads
    )
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.ifEmpty(null))

    //
    // MODULE: mapping
    //
    BWA_MEM(
        CUTADAPT.out.reads,
        PREPARE_GENOME.out.bwa_index
    )
    ch_software_versions = ch_software_versions.mix(BWA_MEM.out.version.ifEmpty(null))

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
    SAMTOOLS_MERGE(mapped_bam)
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_MERGE.out.version.ifEmpty(null))

    //
    // MODULE: mapping stats
    //
    BAM_STAT(SAMTOOLS_MERGE.out.bam)
    ch_software_versions = ch_software_versions.mix(BAM_STAT.out.version.ifEmpty(null))

    //
    // SUBWORKFLOW: filter reads, output pair (like hic pair), raw (pair), and stats
    //
    PAIRTOOLS_PAIRE(
        SAMTOOLS_MERGE.out.bam,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.digest_genome
    )
    ch_software_versions = ch_software_versions.mix(PAIRTOOLS_PAIRE.out.version.ifEmpty(null))

    //
    // combine bin_size and create cooler file, and dump long_bedpe
    //
    cool_bin.combine(PAIRTOOLS_PAIRE.out.pair)
            .map{bin, meta, pair, px -> [meta, bin, pair, px]}
            .set{cool_input}
    COOLER(
        cool_input,
        PREPARE_GENOME.out.chrom_sizes
    )
    ch_software_versions = ch_software_versions.mix(COOLER.out.version.ifEmpty(null))

    //
    // calling ATAC peaks, output ATAC narrowPeak and reads in peak
    //
    ATAC_PEAK(
        PAIRTOOLS_PAIRE.out.raw,
        PREPARE_GENOME.out.chrom_sizes,
        PREPARE_GENOME.out.gsize
    )
    ch_software_versions = ch_software_versions.mix(ATAC_PEAK.out.version.ifEmpty(null))

    //
    // calling distal peaks: [ meta, bin_size, path(macs2), path(long_bedpe), path(short_bed), path(background) ]
    //
    background = MAPS_MULTIENZYME(PREPARE_GENOME.out.fasta, cool_bin, PREPARE_GENOME.out.chrom_sizes).bin_feature
    ch_software_versions = ch_software_versions.mix(MAPS_MULTIENZYME.out.version.ifEmpty(null))
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
    MAPS_PEAK(maps_input)
    ch_software_versions = ch_software_versions.mix(MAPS_PEAK.out.version.ifEmpty(null))

    //
    // Annotate the MAPS peak
    //
    if(!params.skip_peak_annotation){
        MAPS_PEAK.out.peak //[]
            .map{meta, bin_size, peak -> [bin_size, peak]}
            .groupTuple()
            .set{ch_maps_anno}
        BIOC_CHIPPEAKANNO_MAPS(ch_maps_anno, PREPARE_GENOME.out.gtf)
        ch_software_versions = ch_software_versions.mix(BIOC_CHIPPEAKANNO_MAPS.out.version.ifEmpty(null))
        if(!params.skip_virtual_4c){
            BIOC_CHIPPEAKANNO_MAPS.out.csv.mix(COOLER.out.mcool.map{meta, mcool -> [meta.bin, mcool]}.groupTuple())
                                        .groupTuple()
                                        .map{bin, df -> [bin, df[0], df[1]]}
                                        .set{ch_maps_trackviewer}
            //ch_maps_trackviewer.view()
            BIOC_TRACKVIEWER_MAPS(
                ch_maps_trackviewer,
                PAIRTOOLS_PAIRE.out.distalpair.collect{it[1]},
                PREPARE_GENOME.out.gtf,
                PREPARE_GENOME.out.chrom_sizes,
                PREPARE_GENOME.out.digest_genome)
            ch_software_versions = ch_software_versions.mix(BIOC_TRACKVIEWER_MAPS.out.version.ifEmpty(null))
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
            .map{ peak, long_bedpe -> [peak[0], peak[1].flatten(), long_bedpe[1].flatten()] }//bin_size, meta, peak, long_bedpe
            .groupTuple()
            .map{[it[0], it[1].flatten().unique(), it[2].flatten()]}
            .set{ch_diffhicar}
        //ch_diffhicar.view()
        DIFFHICAR(ch_diffhicar)
        ch_software_versions = ch_software_versions.mix(DIFFHICAR.out.version.ifEmpty(null))
        //annotation
        if(!params.skip_peak_annotation){
            BIOC_CHIPPEAKANNO(DIFFHICAR.out.diff, PREPARE_GENOME.out.gtf)
            ch_software_versions = ch_software_versions.mix(BIOC_CHIPPEAKANNO.out.version.ifEmpty(null))
            if(PREPARE_GENOME.out.ucscname) BIOC_ENRICH(BIOC_CHIPPEAKANNO.out.anno.filter{it.size()>0}, PREPARE_GENOME.out.ucscname)
            ch_software_versions = ch_software_versions.mix(BIOC_ENRICH.out.version.ifEmpty(null))
            if(!params.skip_virtual_4c){
                BIOC_CHIPPEAKANNO.out.csv.mix(COOLER.out.mcool.map{meta, mcool -> [meta.bin, mcool]}.groupTuple())
                                            .groupTuple()
                                            .map{bin, df -> [bin, df[0], df[1]]}
                                            .set{ch_trackviewer}
                //ch_trackviewer.view()
                BIOC_TRACKVIEWER(
                    ch_trackviewer,
                    PAIRTOOLS_PAIRE.out.distalpair.collect{it[1]},
                    PREPARE_GENOME.out.gtf,
                    PREPARE_GENOME.out.chrom_sizes,
                    PREPARE_GENOME.out.digest_genome)
                ch_software_versions = ch_software_versions.mix(BIOC_TRACKVIEWER.out.version.ifEmpty(null))
            }
        }
    }

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .flatten()
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )
    if(!params.skip_multiqc){
        //
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowHicar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BAM_STAT.out.idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(PAIRTOOLS_PAIRE.out.stat.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BIOC_CHIPPEAKANNO_MAPS.out.png.collect().ifEmpty([]))
        //ch_multiqc_files = ch_multiqc_files.mix(BIOC_CHIPPEAKANNO.out.png.collect().ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report       = MULTIQC.out.report.toList()
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
