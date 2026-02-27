/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: utils_nfcore_tempo_pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pipeline-specific utility functions for mskcc/tempo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList        } from 'plugin/nf-schema'
include { paramsSummaryLog         } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_INITIALISATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Validate parameters
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required
    //
    if (version) {
        def versionString = ""
        versionString += "mskcc/tempo v${workflow.manifest.version}\n"
        versionString += "Nextflow v${nextflow.version}\n"
        log.info versionString
        System.exit(0)
    }

    //
    // Print parameter summary
    //
    def summary_params = paramsSummaryMap(workflow)

    //
    // Validate input parameters
    //
    if (validate_params) {
        // nf-schema will validate automatically via plugin
    }

    //
    // Check mandatory parameters
    //
    if (!params.input && !params.bamMapping && !(params.aggregate instanceof String && params.aggregate != 'true')) {
        error("Please provide an input samplesheet with --input, a BAM mapping with --bamMapping, or an aggregate TSV with --aggregate")
    }

    //
    // Create channel from input samplesheet
    //
    if (params.input) {
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { row ->
                // nf-schema 2.5.0: row[0] = meta map (fields with "meta" tag),
                // subsequent elements are non-meta fields in schema order
                def meta_raw = row[0]
                def meta = [
                    id:      meta_raw.sample,
                    patient: meta_raw.patient,
                    sample:  meta_raw.sample,
                    status:  meta_raw.status instanceof Integer ? meta_raw.status : meta_raw.status.toInteger(),
                    target:  meta_raw.target,
                    lane:    meta_raw.lane ?: 'L001'
                ]
                def fastq_1 = row[1]
                def fastq_2 = row[2]

                if (!fastq_1) {
                    error("Please check input samplesheet -> FASTQ file for reads 1 does not exist!\n${row}")
                }

                return [ meta, fastq_1, fastq_2 ]
            }
            .set { ch_samplesheet }
    } else {
        ch_samplesheet = Channel.empty()
    }

    //
    // Create channel from BAM mapping (if provided)
    //
    if (params.bamMapping) {
        Channel
            .fromPath(params.bamMapping, checkIfExists: true)
            .splitCsv(sep: '\t', header: true)
            .map { row ->
                if (!row.PATIENT || !row.SAMPLE || !row.TARGET || !row.BAM || !row.BAI) {
                    error("bamMapping TSV must have columns: PATIENT, SAMPLE, TARGET, BAM, BAI (STATUS is optional, defaults to 0). Found: ${row.keySet()}")
                }
                def meta = [
                    id:      row.SAMPLE,
                    patient: row.PATIENT,
                    sample:  row.SAMPLE,
                    status:  row.STATUS ? row.STATUS.toInteger() : 0,
                    target:  row.TARGET
                ]
                def bam = file(row.BAM, checkIfExists: true)
                def bai = file(row.BAI, checkIfExists: true)
                if (!bam.name.endsWith('.bam')) {
                    error("BAM file must end with .bam: ${row.BAM}")
                }
                if (!bai.name.endsWith('.bai') && !bai.name.endsWith('.bam.bai')) {
                    error("BAI file must end with .bai: ${row.BAI}")
                }
                return [ meta, bam, bai ]
            }
            .set { ch_bam_input }
    } else {
        ch_bam_input = Channel.empty()
    }

    emit:
    samplesheet = ch_samplesheet
    bam_input   = ch_bam_input
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: PIPELINE_COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    versions        // channel: All software versions

    main:
    // Completion handlers are registered in the entry workflow scope (main.nf)
    log.debug "Pipeline completion subworkflow invoked"
}
