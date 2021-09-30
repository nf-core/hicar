#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/hicar
========================================================================================
    Github : https://github.com/nf-core/hicar
    Website: https://nf-co.re/hicar
    Slack  : https://nfcore.slack.com/channels/hicar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta      = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.bwa_index  = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.gtf        = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff        = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed   = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.macs_gsize = WorkflowMain.getGenomeAttribute(params, 'macs_gsize')
anno_readme       = WorkflowMain.getGenomeAttribute(params, 'readme')
// Save AWS IGenomes file containing annotation version
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { HICAR } from './workflows/hicar'

//
// WORKFLOW: Run main nf-core/hicar analysis pipeline
//
workflow NFCORE_HICAR {
    HICAR ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_HICAR ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
