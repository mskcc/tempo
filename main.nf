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
include { BWAMEM2_INDEX           } from './modules/local/bwamem2/index/main'
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
    bam_input   // channel: BAM inputs from --bamMapping

    main:

    //
    // WORKFLOW: Run pipeline
    //

    // ---------------------------------------------------------------
    // Resolve genome references at RUNTIME (not config-parse time).
    //
    // CLI params (--genome, --reference_base) are only guaranteed
    // available here, not during config parsing.  If flat nf-core
    // params (e.g. params.fasta) are already set (via -profile test
    // or CLI --fasta), they take priority.
    //
    // NOTE: params is READ-ONLY in Nextflow 24+, so we compute into
    // local variables (ref_*) and use those for channel creation.
    // ---------------------------------------------------------------
    def ref_fasta             = params.fasta
    def ref_fasta_fai         = params.fasta_fai
    def ref_dict              = params.dict
    def ref_bwa_index         = params.bwa_index
    def ref_dbsnp             = params.dbsnp
    def ref_dbsnp_tbi         = params.dbsnp_tbi
    def ref_known_indels      = params.known_indels
    def ref_known_indels_tbi  = params.known_indels_tbi
    def ref_germline_resource     = params.germline_resource
    def ref_germline_resource_tbi = params.germline_resource_tbi
    def ref_intervals         = params.intervals
    def ref_pon               = params.pon
    def ref_pon_tbi           = params.pon_tbi

    if (params.genome && params.reference_base && !ref_fasta) {
        def rb = params.reference_base.replaceAll('/+$', '')  // strip trailing slash
        def genome_base = params.genome == 'GRCh37'
            ? "${rb}/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37"
            : params.genome == 'GRCh38'
                ? "${rb}/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh38"
                : "${rb}/mskcc-igenomes/igenomes/smallGRCh37"
        def targets_base = "${rb}/mskcc-igenomes/${params.genome.toLowerCase()}/tempo_targets"

        // Core GATK references
        ref_fasta             = "${genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
        ref_fasta_fai         = "${ref_fasta}.fai"
        ref_dict              = "${genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict"
        ref_bwa_index         = "${genome_base}/Sequence/BWAIndex/human_g1k_v37_decoy.fasta"
        ref_dbsnp             = "${genome_base}/Annotation/GATKBundle/dbsnp_138.b37.vcf"
        ref_dbsnp_tbi         = "${ref_dbsnp}.idx"
        ref_known_indels      = "${genome_base}/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf"
        ref_known_indels_tbi  = "${genome_base}/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf.idx"
        ref_intervals         = ref_intervals ?: "${genome_base}/Annotation/intervals/human.b37.genome.bed"

        // Tempo-specific references
        ref_germline_resource     = ref_germline_resource     ?: "${rb}/mskcc-igenomes/grch37/gnomad/gnomad.exomes.r2.1.1.sites.non_cancer.vcf.gz"
        ref_germline_resource_tbi = ref_germline_resource_tbi ?: "${ref_germline_resource}.tbi"

        // PON: depends on assay type
        if (!ref_pon) {
            if (params.assay_type == 'genome') {
                ref_pon     = "${rb}/mskcc-igenomes/grch37/annotation/wgs.pon.vcf.gz"
                ref_pon_tbi = "${ref_pon}.tbi"
            } else {
                ref_pon     = "${rb}/mskcc-igenomes/grch37/annotation/wes.pon.vcf.gz"
                ref_pon_tbi = "${ref_pon}.tbi"
            }
        }

        log.info "Resolved references for genome '${params.genome}' from reference_base: ${rb}"
    }

    // Validate required reference parameters — fail early with clear message
    def ref_map = [
        fasta: ref_fasta, fasta_fai: ref_fasta_fai, dict: ref_dict,
        bwa_index: ref_bwa_index, dbsnp: ref_dbsnp, dbsnp_tbi: ref_dbsnp_tbi,
        known_indels: ref_known_indels, known_indels_tbi: ref_known_indels_tbi
    ]
    def missing = ref_map.findAll { k, v -> !v }.collect { k, v -> k }
    if (missing) {
        error "Missing required reference parameter(s): ${missing.join(', ')}. " +
              "Please provide --genome and --reference_base, use -profile test, " +
              "or set each reference path individually (--fasta, --dbsnp, etc.)."
    }

    // Prepare reference channels as value channels (nf-core/sarek pattern)
    // Using Channel.value() ensures these pair correctly with every sample in multi-sample runs
    ch_fasta            = Channel.value([ [id:'genome'], file(ref_fasta, checkIfExists: true) ])
    ch_fasta_fai        = Channel.value([ [id:'genome'], file(ref_fasta_fai, checkIfExists: true) ])
    ch_dict             = Channel.value([ [id:'genome'], file(ref_dict, checkIfExists: true) ])
    // bwa-mem2 index: check if pre-built index exists, otherwise build on-the-fly
    def bwamem2_index_exists = ref_bwa_index ? file("${ref_bwa_index}.bwt.2bit.64").exists() : false
    if (ref_bwa_index && bwamem2_index_exists) {
        // Pre-built bwa-mem2 index available
        ch_bwa_index = Channel.fromPath("${ref_bwa_index}.{amb,ann,bwt.2bit.64,pac,0123}", checkIfExists: true)
            .collect()
            .map { files -> [ [id:'genome'], files ] }
    } else if (ref_bwa_index) {
        // No bwa-mem2 index — build on-the-fly from fasta
        BWAMEM2_INDEX ( ch_fasta )
        ch_bwa_index = BWAMEM2_INDEX.out.index
    }
    ch_dbsnp            = Channel.value([ [id:'dbsnp'], file(ref_dbsnp, checkIfExists: true) ])
    ch_dbsnp_tbi        = Channel.value([ [id:'dbsnp'], file(ref_dbsnp_tbi, checkIfExists: true) ])
    ch_known_indels     = Channel.value([ [id:'indels'], file(ref_known_indels, checkIfExists: true) ])
    ch_known_indels_tbi = Channel.value([ [id:'indels'], file(ref_known_indels_tbi, checkIfExists: true) ])
    ch_germline_resource     = ref_germline_resource     ? Channel.value([ [id:'gnomad'], file(ref_germline_resource, checkIfExists: true) ])     : Channel.value([ [id:'gnomad'], [] ])
    ch_germline_resource_tbi = ref_germline_resource_tbi ? Channel.value([ [id:'gnomad'], file(ref_germline_resource_tbi, checkIfExists: true) ]) : Channel.value([ [id:'gnomad'], [] ])
    ch_intervals        = ref_intervals        ? Channel.value([ file(ref_intervals, checkIfExists: true) ])                       : Channel.value([])
    ch_pon              = ref_pon              ? Channel.value([ [id:'pon'], file(ref_pon, checkIfExists: true) ])                  : Channel.value([ [id:'pon'], [] ])
    ch_pon_tbi          = ref_pon_tbi          ? Channel.value([ [id:'pon'], file(ref_pon_tbi, checkIfExists: true) ])              : Channel.value([ [id:'pon'], [] ])

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
        ch_pon_tbi,
        bam_input
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
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.bam_input
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
