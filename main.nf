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
    // ---------------------------------------------------------------
    if (params.genome && params.reference_base && !params.fasta) {
        def rb = params.reference_base
        def genome_base = params.genome == 'GRCh37'
            ? "${rb}/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37"
            : params.genome == 'GRCh38'
                ? "${rb}/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh38"
                : "${rb}/mskcc-igenomes/igenomes/smallGRCh37"
        def targets_base = "${rb}/mskcc-igenomes/${params.genome.toLowerCase()}/tempo_targets"

        // Core GATK references
        params.fasta             = "${genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
        params.fasta_fai         = "${params.fasta}.fai"
        params.dict              = "${genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict"
        params.bwa_index         = "${genome_base}/Sequence/BWAIndex/human_g1k_v37_decoy.fasta.{amb,ann,bwt,pac,sa}"
        params.dbsnp             = "${genome_base}/Annotation/GATKBundle/dbsnp_138.b37.vcf"
        params.dbsnp_tbi         = "${params.dbsnp}.idx"
        params.known_indels      = "${genome_base}/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf"
        params.known_indels_tbi  = "${genome_base}/Annotation/GATKBundle/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.b37.vcf.idx"
        params.intervals         = params.intervals ?: "${genome_base}/Annotation/intervals/human.b37.genome.bed"

        // Tempo-specific references
        params.germline_resource     = params.germline_resource     ?: "${rb}/mskcc-igenomes/grch37/gnomad/gnomad.exomes.r2.1.1.sites.non_cancer.vcf.gz"
        params.germline_resource_tbi = params.germline_resource_tbi ?: "${params.germline_resource}.tbi"
        params.msi_sensor_list       = params.msi_sensor_list       ?: "${genome_base}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.microsatellites.list"
        params.facets_vcf            = params.facets_vcf            ?: "${rb}/mskcc-igenomes/igenomes/Homo_sapiens/GATK/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf"
        params.delly_exclude_regions = params.delly_exclude_regions ?: "${rb}/mskcc-igenomes/grch37/delly/human.hg19.excl.tsv"
        params.snp_gc_corrections    = params.snp_gc_corrections    ?: "${rb}/mskcc-igenomes/grch37/ascat/SnpGcCorrections.tsv"

        // PON: depends on assay type
        if (!params.pon) {
            if (params.assay_type == 'genome') {
                params.pon     = "${rb}/mskcc-igenomes/grch37/annotation/wgs.pon.vcf.gz"
                params.pon_tbi = "${params.pon}.tbi"
            } else {
                params.pon     = "${rb}/mskcc-igenomes/grch37/annotation/wes.pon.vcf.gz"
                params.pon_tbi = "${params.pon}.tbi"
            }
        }

        // Additional Tempo references
        params.bait_intervals    = params.bait_intervals    ?: "${targets_base}/\${targets_id}/baits.interval_list"
        params.target_intervals  = params.target_intervals  ?: "${targets_base}/\${targets_id}/targets.interval_list"
        params.splice_sites      = params.splice_sites      ?: "${rb}/mskcc-igenomes/grch37/splice_sites/splice_sites.bed"
        params.hla_fasta         = params.hla_fasta         ?: "${rb}/mskcc-igenomes/grch37/hla/abc_complete.fasta"
        params.hla_dat           = params.hla_dat           ?: "${rb}/mskcc-igenomes/grch37/hla/hla.dat"
        params.neoantigen_cdna   = params.neoantigen_cdna   ?: "${rb}/mskcc-igenomes/grch37/neoantigen/Homo_sapiens.GRCh37.75.cdna.all.fa.gz"
        params.neoantigen_cds    = params.neoantigen_cds    ?: "${rb}/mskcc-igenomes/grch37/neoantigen/Homo_sapiens.GRCh37.75.cds.all.fa.gz"
        params.vep_cache         = params.vep_cache         ?: "${rb}/mskcc-igenomes/grch37/vep"

        log.info "Resolved references for genome '${params.genome}' from reference_base: ${rb}"
    }

    // Validate required reference parameters — fail early with clear message
    def required_refs = [
        'fasta', 'fasta_fai', 'dict', 'bwa_index',
        'dbsnp', 'dbsnp_tbi', 'known_indels', 'known_indels_tbi'
    ]
    def missing = required_refs.findAll { !params[it] }
    if (missing) {
        error "Missing required reference parameter(s): ${missing.join(', ')}. " +
              "Please provide --genome and --reference_base, use -profile test, " +
              "or set each reference path individually (--fasta, --dbsnp, etc.)."
    }

    // Prepare reference channels as value channels (nf-core/sarek pattern)
    // Using Channel.value() ensures these pair correctly with every sample in multi-sample runs
    ch_fasta            = Channel.value([ [id:'genome'], file(params.fasta, checkIfExists: true) ])
    ch_fasta_fai        = Channel.value([ [id:'genome'], file(params.fasta_fai, checkIfExists: true) ])
    ch_dict             = Channel.value([ [id:'genome'], file(params.dict, checkIfExists: true) ])
    // bwa-mem2 index: check if pre-built index exists, otherwise build on-the-fly
    def bwamem2_index_exists = params.bwa_index ? file("${params.bwa_index}.bwt.2bit.64").exists() : false
    if (params.bwa_index && bwamem2_index_exists) {
        // Pre-built bwa-mem2 index available
        ch_bwa_index = Channel.fromPath("${params.bwa_index}.{amb,ann,bwt.2bit.64,pac,0123}", checkIfExists: true)
            .collect()
            .map { files -> [ [id:'genome'], files ] }
    } else if (params.bwa_index) {
        // No bwa-mem2 index — build on-the-fly from fasta
        BWAMEM2_INDEX ( ch_fasta )
        ch_bwa_index = BWAMEM2_INDEX.out.index
    }
    ch_dbsnp            = Channel.value([ [id:'dbsnp'], file(params.dbsnp, checkIfExists: true) ])
    ch_dbsnp_tbi        = Channel.value([ [id:'dbsnp'], file(params.dbsnp_tbi, checkIfExists: true) ])
    ch_known_indels     = Channel.value([ [id:'indels'], file(params.known_indels, checkIfExists: true) ])
    ch_known_indels_tbi = Channel.value([ [id:'indels'], file(params.known_indels_tbi, checkIfExists: true) ])
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
