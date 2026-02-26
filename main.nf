#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mskcc/tempo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mskcc/tempo
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TEMPO                   } from './workflows/tempo'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_tempo_pipeline/main'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_tempo_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow MSKCC_TEMPO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //

    // Prepare reference channels as value channels (nf-core/sarek pattern)
    // Using Channel.value() ensures these pair correctly with every sample in multi-sample runs
    ch_fasta            = params.fasta            ? Channel.value([ [id:'genome'], file(params.fasta, checkIfExists: true) ])            : Channel.value([ [id:'genome'], [] ])
    ch_fasta_fai        = params.fasta_fai        ? Channel.value([ [id:'genome'], file(params.fasta_fai, checkIfExists: true) ])        : Channel.value([ [id:'genome'], [] ])
    ch_dict             = params.dict             ? Channel.value([ [id:'genome'], file(params.dict, checkIfExists: true) ])             : Channel.value([ [id:'genome'], [] ])
    ch_bwa_index        = params.bwa_index        ? Channel.fromPath("${params.bwa_index}.{amb,ann,bwt,pac,sa,0123,bwt.2bit.64}", checkIfExists: false).collect().map{ files -> [ [id:'genome'], files ] }        : Channel.value([ [id:'genome'], [] ])
    ch_dbsnp            = params.dbsnp            ? Channel.value([ [id:'dbsnp'], file(params.dbsnp, checkIfExists: true) ])             : Channel.value([ [id:'dbsnp'], [] ])
    ch_dbsnp_tbi        = params.dbsnp_tbi        ? Channel.value([ [id:'dbsnp'], file(params.dbsnp_tbi, checkIfExists: true) ])         : Channel.value([ [id:'dbsnp'], [] ])
    ch_known_indels     = params.known_indels     ? Channel.value([ [id:'indels'], file(params.known_indels, checkIfExists: true) ])     : Channel.value([ [id:'indels'], [] ])
    ch_known_indels_tbi = params.known_indels_tbi ? Channel.value([ [id:'indels'], file(params.known_indels_tbi, checkIfExists: true) ]) : Channel.value([ [id:'indels'], [] ])
    ch_germline_resource     = params.germline_resource     ? Channel.value([ [id:'gnomad'], file(params.germline_resource, checkIfExists: true) ])     : Channel.value([ [id:'gnomad'], [] ])
    ch_germline_resource_tbi = params.germline_resource_tbi ? Channel.value([ [id:'gnomad'], file(params.germline_resource_tbi, checkIfExists: true) ]) : Channel.value([ [id:'gnomad'], [] ])
    ch_intervals        = params.intervals        ? Channel.value([ file(params.intervals, checkIfExists: true) ])                       : Channel.value([])
    ch_pon              = params.pon              ? Channel.value([ [id:'pon'], file(params.pon, checkIfExists: true) ])                  : Channel.value([ [id:'pon'], [] ])
    ch_pon_tbi          = params.pon_tbi          ? Channel.value([ [id:'pon'], file(params.pon_tbi, checkIfExists: true) ])              : Channel.value([ [id:'pon'], [] ])

    TEMPO (
        samplesheet,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_bwa_index,
        ch_dbsnp,
        ch_dbsnp_tbi,
        ch_known_indels,
        ch_known_indels_tbi,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_intervals,
        ch_pon,
        ch_pon_tbi
    )

    emit:
    versions = TEMPO.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Register completion handlers at global scope (outside workflow blocks)
workflow.onComplete {
    if (workflow.success) {
        log.info "-[mskcc/tempo] Pipeline completed successfully-"
    } else {
        log.info "-[mskcc/tempo] Pipeline completed with errors-"
    }
}

workflow.onError {
    log.error "Pipeline failed. See .nextflow.log for details."
}

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir
    )

    //
    // WORKFLOW: Run main workflow
    //
    MSKCC_TEMPO (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        MSKCC_TEMPO.out.versions
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
