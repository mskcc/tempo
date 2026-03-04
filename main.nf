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
    // Resolve genome references: bridge original Tempo params.genomes
    // structure (from conf/references.config) to nf-core flat params.
    // If flat params (e.g. params.fasta) are already set (e.g. via
    // -profile test), they take priority.  Otherwise, look up from
    // params.genomes[params.genome].
    // ---------------------------------------------------------------
    def genomeRef = (params.genome && params.genomes && params.genomes.containsKey(params.genome))
                    ? params.genomes[params.genome]
                    : [:]

    // Helper: resolve a flat param from the genomes map, with a mapping
    // from nf-core flat param name -> original Tempo genomes key name
    def refMapping = [
        fasta            : 'genomeFile',
        fasta_fai        : 'genomeIndex',
        dict             : 'genomeDict',
        bwa_index        : 'bwaIndex',
        dbsnp            : 'dbsnp',
        dbsnp_tbi        : 'dbsnpIndex',
        known_indels     : 'knownIndels',
        known_indels_tbi : 'knownIndelsIndex',
        intervals        : 'intervals',
        germline_resource     : 'gnomadWesVcf',
        germline_resource_tbi : 'gnomadWesVcfIndex',
        msi_sensor_list  : 'msiSensorList',
        facets_vcf       : 'facetsVcf',
        delly_exclude_regions : 'svCallingExcludeRegions',
        snp_gc_corrections    : 'snpGcCorrections',
    ]

    // For each mapping, if the flat param is not set, fill from genomes map
    refMapping.each { flatKey, genomesKey ->
        if (!params[flatKey] && genomeRef[genomesKey]) {
            params[flatKey] = genomeRef[genomesKey]
        }
    }

    // Resolve PON based on assay_type (exome vs genome)
    if (!params.pon && genomeRef) {
        def ponKey = params.assay_type == 'genome' ? 'wgsPoN' : 'exomePoN'
        def ponIdxKey = params.assay_type == 'genome' ? 'wgsPoNIndex' : 'exomePoNIndex'
        if (genomeRef[ponKey])    params.pon     = genomeRef[ponKey]
        if (genomeRef[ponIdxKey]) params.pon_tbi = genomeRef[ponIdxKey]
    }

    // Validate required reference parameters — fail early with clear message
    def required_refs = [
        'fasta', 'fasta_fai', 'dict', 'bwa_index',
        'dbsnp', 'dbsnp_tbi', 'known_indels', 'known_indels_tbi'
    ]
    def missing = required_refs.findAll { !params[it] }
    if (missing) {
        error "Missing required reference parameter(s): ${missing.join(', ')}. " +
              "Please provide all reference files via params, --genome with conf/references.config, " +
              "or a config profile (e.g., -profile test)."
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
