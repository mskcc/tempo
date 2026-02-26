/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed from nf-core/modules
//
include { FASTQC                                       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                      } from '../modules/nf-core/multiqc/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { BWAMEM2_MEM                                  } from '../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_SORT                                } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SORTED       } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MD           } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL        } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE                               } from '../modules/nf-core/samtools/merge/main'
include { GATK4_MARKDUPLICATES                         } from '../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR                       } from '../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR                              } from '../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_MUTECT2                                } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR  } from '../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL } from '../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION                 } from '../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_LEARNREADORIENTATIONMODEL              } from '../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_FILTERMUTECTCALLS                      } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_HAPLOTYPECALLER                        } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { STRELKA_SOMATIC                              } from '../modules/nf-core/strelka/somatic/main'
include { MANTA_SOMATIC                                } from '../modules/nf-core/manta/somatic/main'
include { ENSEMBLVEP_VEP                               } from '../modules/nf-core/ensemblvep/vep/main'
include { MSISENSORPRO_SCAN                            } from '../modules/nf-core/msisensorpro/scan/main'
include { MSISENSORPRO_MSISOMATIC                      } from '../modules/nf-core/msisensorpro/msisomatic/main'
include { PICARD_COLLECTHSMETRICS                      } from '../modules/nf-core/picard/collecthsmetrics/main'
include { QUALIMAP_BAMQC                               } from '../modules/nf-core/qualimap/bamqc/main'

//
// MODULE: Local modules (Tempo-specific, no nf-core equivalent)
//
include { POLYSOLVER                  } from '../modules/local/polysolver/main'
include { LOHHLA                      } from '../modules/local/lohhla/main'
include { SNPPILEUP                   } from '../modules/local/snppileup/main'
include { FACETS                      } from '../modules/local/facets/main'
include { CONPAIR_PILEUP              } from '../modules/local/conpair/pileup/main'
include { CONPAIR_CONCORDANCE         } from '../modules/local/conpair/concordance/main'
include { VCF2MAF                     } from '../modules/local/vcf2maf/main'
include { DELLY_CALL_SOMATIC          } from '../modules/local/delly/call/main'
include { DELLY_COMBINE               } from '../modules/local/delly/combine/main'
include { DELLY_CALL_GERMLINE         } from '../modules/local/delly/call_germline/main'
include { STRELKA2_COMBINE_SOMATIC    } from '../modules/local/strelka2/combine_somatic/main'
include { SOMATIC_COMBINE_CHANNEL     } from '../modules/local/somatic/combine_channel/main'
include { SOMATIC_MERGE_SV            } from '../modules/local/somatic/merge_sv/main'
include { STRELKA2_GERMLINE           } from '../modules/local/strelka2/germline/main'
include { MANTA_GERMLINE              } from '../modules/local/manta/germline/main'
include { GERMLINE_MERGE_SV           } from '../modules/local/germline/merge_sv/main'
include { GERMLINE_COMBINE_CHANNEL    } from '../modules/local/germline/combine_channel/main'

// New Phase 2 local modules — SV pipeline
include { SVABA_SOMATIC               } from '../modules/local/svaba/somatic/main'
include { SVABA_GERMLINE              } from '../modules/local/svaba/germline/main'
include { SVTOOLS_VCF2BEDPE_SOMATIC   } from '../modules/local/svtools/vcf2bedpe_somatic/main'
include { SVTOOLS_VCF2BEDPE_GERMLINE  } from '../modules/local/svtools/vcf2bedpe_germline/main'
include { IANNOTATESV_SOMATIC         } from '../modules/local/iannotatesv/somatic/main'
include { IANNOTATESV_GERMLINE        } from '../modules/local/iannotatesv/germline/main'
include { CLUSTERSV                   } from '../modules/local/clustersv/main'
include { SVCIRCOS                    } from '../modules/local/svcircos/main'
include { SVCLONE                     } from '../modules/local/svclone/main'

// New Phase 2 local modules — WGS-only (ASCAT + BRASS + HRDetect)
include { ASCAT_ALLELECOUNT           } from '../modules/local/ascat/allelecount/main'
include { ASCAT_RUN                   } from '../modules/local/ascat/run/main'
include { BRASS_GENERATE_BAS          } from '../modules/local/brass/generate_bas/main'
include { BRASS_INPUT                 } from '../modules/local/brass/input/main'
include { BRASS_COVER                 } from '../modules/local/brass/cover/main'
include { BRASS_RUN                   } from '../modules/local/brass/run/main'
include { HRDETECT                    } from '../modules/local/hrdetect/main'
include { SV_SIGNATURES               } from '../modules/local/svsignatures/main'

// New Phase 2 local modules — Annotation & Signatures
include { GERMLINE_ANNOTATE_MAF       } from '../modules/local/germline/annotate_maf/main'
include { SOMATIC_FACETS_ANNOTATION   } from '../modules/local/somatic/facets_annotation/main'
include { GERMLINE_FACETS_ANNOTATION  } from '../modules/local/germline/facets_annotation/main'
include { FACETS_PREVIEW_QC           } from '../modules/local/facets/preview_qc/main'
include { NEOANTIGEN                  } from '../modules/local/neoantigen/main'
include { MUTSIG                      } from '../modules/local/mutsig/main'
include { GERMLINE_COMBINE_HC_VCF     } from '../modules/local/germline/combine_hc_vcf/main'
include { SPLIT_INTERVALS             } from '../modules/local/splitintervals/main'
include { METADATA_PARSER             } from '../modules/local/metadata_parser/main'

// New Phase 2 local modules — QC & Reporting
include { ALFRED                      } from '../modules/local/alfred/main'
include { CONPAIR_ALL                 } from '../modules/local/conpair/all/main'
include { MULTIQC_SAMPLE              } from '../modules/local/multiqc/sample/main'
include { MULTIQC_SOMATIC             } from '../modules/local/multiqc/somatic/main'
include { MULTIQC_COHORT              } from '../modules/local/multiqc/cohort/main'

// New Phase 2 local modules — Aggregation
include { AGGREGATE_SOMATIC_MAF         } from '../modules/local/aggregate/somatic_maf/main'
include { AGGREGATE_SOMATIC_SV          } from '../modules/local/aggregate/somatic_sv/main'
include { AGGREGATE_SOMATIC_FACETS      } from '../modules/local/aggregate/somatic_facets/main'
include { AGGREGATE_SOMATIC_NETMHC      } from '../modules/local/aggregate/somatic_netmhc/main'
include { AGGREGATE_SOMATIC_METADATA    } from '../modules/local/aggregate/somatic_metadata/main'
include { AGGREGATE_SOMATIC_LOHHLA      } from '../modules/local/aggregate/somatic_lohhla/main'
include { AGGREGATE_SOMATIC_HRDETECT    } from '../modules/local/aggregate/somatic_hrdetect/main'
include { AGGREGATE_SOMATIC_SVCLONE     } from '../modules/local/aggregate/somatic_svclone/main'
include { AGGREGATE_SOMATIC_SVSIGNATURES } from '../modules/local/aggregate/somatic_svsignatures/main'
include { AGGREGATE_GERMLINE_MAF        } from '../modules/local/aggregate/germline_maf/main'
include { AGGREGATE_GERMLINE_SV         } from '../modules/local/aggregate/germline_sv/main'
include { AGGREGATE_QC_BAM             } from '../modules/local/aggregate/qc_bam/main'
include { AGGREGATE_QC_CONPAIR         } from '../modules/local/aggregate/qc_conpair/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TEMPO {

    take:
    ch_input            // channel: samplesheet parsed rows [ meta, fastq_1, fastq_2 ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fasta_fai        // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    ch_bwa_index        // channel: [ val(meta), path(index) ]
    ch_dbsnp            // channel: [ val(meta), path(vcf) ]
    ch_dbsnp_tbi        // channel: [ val(meta), path(tbi) ]
    ch_known_indels     // channel: [ val(meta), path(vcf) ]
    ch_known_indels_tbi // channel: [ val(meta), path(tbi) ]
    ch_germline_resource     // channel: [ val(meta), path(vcf) ]
    ch_germline_resource_tbi // channel: [ val(meta), path(tbi) ]
    ch_intervals        // channel: [ path(intervals) ]
    ch_pon              // channel: [ val(meta), path(vcf) ]
    ch_pon_tbi          // channel: [ val(meta), path(tbi) ]

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // =============================================
    // WORKFLOW CONTROL FLAGS (mirrors original Tempo dsl2.nf)
    // =============================================

    // Parse --workflows string into boolean flags
    def WFs = params.workflows instanceof Boolean ? '' : (params.workflows ?: '')
    def wfList = WFs.split(',').collect{ it.trim().toLowerCase() }.unique().findAll{ it }

    def doWF_SNV        = ['snv', 'mutsig'].any{ it in wfList }
    def doWF_SV         = 'sv' in wfList
    def doWF_manta      = ['snv', 'sv', 'mutsig'].any{ it in wfList }
    def doWF_facets     = ['lohhla', 'facets', 'snv', 'mutsig', 'germsnv'].any{ it in wfList }
    def doWF_loh        = ['lohhla', 'snv', 'mutsig'].any{ it in wfList }
    def doWF_germSNV    = 'germsnv' in wfList
    def doWF_germSV     = 'germsv' in wfList
    def doWF_QC         = 'qc' in wfList
    def doWF_msiSensor  = 'msisensor' in wfList
    def doWF_mutSig     = 'mutsig' in wfList
    def doWF_mdParse    = doWF_manta && doWF_facets && doWF_loh && doWF_SNV && doWF_msiSensor && doWF_mutSig

    // WGS-only features gated by assay_type
    def isWGS = params.assay_type == 'genome'

    // If SV workflow + WGS + facets-based CNV source, ensure facets runs
    if (doWF_SV && isWGS && ['hisens','purity'].contains(params.svcnv)) {
        doWF_facets = true
    }

    //
    // Parse samplesheet and group by patient/sample
    //
    ch_input
        .map { meta, fastq_1, fastq_2 ->
            def new_meta = meta + [
                id:         meta.sample,
                read_group: "@RG\\tID:${meta.sample}_${meta.lane}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}"
            ]
            [ new_meta, [ fastq_1, fastq_2 ] ]
        }
        .set { ch_reads }

    // =============================================
    // RAW READ QC
    // =============================================

    //
    // MODULE: FastQC
    // input: tuple val(meta), path(reads)
    //
    FASTQC ( ch_reads )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    // =============================================
    // PREPROCESSING: TRIM + ALIGN + SORT
    // =============================================

    //
    // MODULE: fastp
    // input: tuple val(meta), path(reads), path(adapter_fasta)
    //        val discard_trimmed_pass
    //        val save_trimmed_fail
    //        val save_merged
    //
    ch_reads
        .map { meta, reads -> [ meta, reads, [] ] }
        .set { ch_fastp_input }

    FASTP (
        ch_fastp_input,
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
        false   // save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})

    //
    // MODULE: BWA-MEM2
    // input: tuple val(meta), path(reads)
    //        tuple val(meta2), path(index)
    //        tuple val(meta3), path(fasta)
    //        val sort_bam
    //
    BWAMEM2_MEM (
        FASTP.out.reads,
        ch_bwa_index,
        ch_fasta,
        true   // sort_bam
    )

    //
    // MODULE: Samtools sort
    // input: tuple val(meta), path(bam)
    //        tuple val(meta2), path(fasta)
    //        val index_format
    //
    SAMTOOLS_SORT (
        BWAMEM2_MEM.out.bam,
        ch_fasta,
        []  // index_format - use default
    )

    //
    // MODULE: Samtools index
    // input: tuple val(meta), path(input)
    //
    SAMTOOLS_INDEX_SORTED ( SAMTOOLS_SORT.out.bam )

    // =============================================
    // MULTI-LANE MERGE
    // =============================================

    //
    // Group BAMs by sample for merging (multi-lane)
    //
    SAMTOOLS_SORT.out.bam
        .map { meta, bam ->
            def new_meta = meta.subMap('patient', 'sample', 'status', 'sex') + [id: meta.sample]
            [ new_meta, bam ]
        }
        .groupTuple()
        .branch {
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }
        .set { ch_bams_to_merge }

    //
    // MODULE: Samtools merge
    // input: tuple val(meta), path(input_files, stageAs: "?/*")
    //        tuple val(meta2), path(fasta), path(fai), path(gzi)
    //
    ch_fasta
        .combine(ch_fasta_fai.map{ it[1] })
        .map { meta, fasta, fai -> [ meta, fasta, fai, [] ] }
        .first()
        .set { ch_ref_fasta_fai_gzi }

    SAMTOOLS_MERGE (
        ch_bams_to_merge.multiple,
        ch_ref_fasta_fai_gzi
    )

    // Combine single-lane and merged BAMs
    ch_bams_to_merge.single
        .map { meta, bam -> [ meta, bam[0] ] }
        .mix(SAMTOOLS_MERGE.out.bam)
        .set { ch_bams_merged }

    // =============================================
    // MARK DUPLICATES
    // =============================================

    //
    // MODULE: GATK4 MarkDuplicates
    // input: tuple val(meta), path(bam)
    //        path fasta
    //        path fasta_fai
    //
    GATK4_MARKDUPLICATES (
        ch_bams_merged,
        ch_fasta.map{ it[1] },
        ch_fasta_fai.map{ it[1] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})

    // Index the marked BAMs
    SAMTOOLS_INDEX_MD ( GATK4_MARKDUPLICATES.out.bam )

    // =============================================
    // BASE QUALITY SCORE RECALIBRATION
    // =============================================

    //
    // MODULE: GATK4 BaseRecalibrator
    // input: tuple val(meta), path(input), path(input_index), path(intervals)
    //        tuple val(meta2), path(fasta)
    //        tuple val(meta3), path(fai)
    //        tuple val(meta4), path(dict)
    //        tuple val(meta5), path(known_sites)
    //        tuple val(meta6), path(known_sites_tbi)
    //
    GATK4_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX_MD.out.bai)
        .map { meta, bam, bai -> [ meta, bam, bai, [] ] }  // empty intervals
        .set { ch_md_bam_bai_intervals }

    // Combine known sites: dbsnp + known_indels
    ch_dbsnp
        .map { meta, vcf -> vcf }
        .mix(ch_known_indels.map { meta, vcf -> vcf })
        .collect()
        .map { files -> [ [id: 'known_sites'], files ] }
        .set { ch_known_sites }

    ch_dbsnp_tbi
        .map { meta, tbi -> tbi }
        .mix(ch_known_indels_tbi.map { meta, tbi -> tbi })
        .collect()
        .map { files -> [ [id: 'known_sites_tbi'], files ] }
        .set { ch_known_sites_tbi }

    GATK4_BASERECALIBRATOR (
        ch_md_bam_bai_intervals,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )

    //
    // MODULE: GATK4 ApplyBQSR
    // input: tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    //        path fasta
    //        path fai
    //        path dict
    //
    GATK4_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX_MD.out.bai)
        .join(GATK4_BASERECALIBRATOR.out.table)
        .map { meta, bam, bai, table -> [ meta, bam, bai, table, [] ] }  // empty intervals
        .set { ch_bam_bai_table_intervals }

    GATK4_APPLYBQSR (
        ch_bam_bai_table_intervals,
        ch_fasta.map{ it[1] },
        ch_fasta_fai.map{ it[1] },
        ch_dict.map{ it[1] }
    )

    // Index recalibrated BAMs
    SAMTOOLS_INDEX_RECAL ( GATK4_APPLYBQSR.out.bam )

    // Create channel of recalibrated BAMs with index
    GATK4_APPLYBQSR.out.bam
        .join(SAMTOOLS_INDEX_RECAL.out.bai)
        .set { ch_recal_bam_bai }

    // =============================================
    // TUMOR-NORMAL PAIRING
    // =============================================

    ch_recal_bam_bai
        .branch {
            tumor:  it[0].status == 1
            normal: it[0].status == 0
        }
        .set { ch_recal_branched }

    // Create tumor-normal pairs by patient
    ch_recal_branched.tumor
        .map { meta, bam, bai -> [ meta.patient, meta, bam, bai ] }
        .combine(
            ch_recal_branched.normal.map { meta, bam, bai -> [ meta.patient, meta, bam, bai ] },
            by: 0
        )
        .map { patient, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            def pair_meta = [
                id:        "${tumor_meta.sample}__${normal_meta.sample}",
                patient:   patient,
                tumor_id:  tumor_meta.sample,
                normal_id: normal_meta.sample,
                status:    1,
                sex:       tumor_meta.sex ?: 'NA'
            ]
            [ pair_meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
        }
        .set { ch_tumor_normal_pair }

    // =============================================
    // SOMATIC SNV/INDEL CALLING
    // =============================================

    // =============================================
    // MANTA SOMATIC (needed by both SNV and SV workflows)
    // =============================================

    if (doWF_manta) {
        //
        // MODULE: Manta somatic (provides indel candidates for Strelka2 + SV calls)
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, nbam, nbai, tbam, tbai, [], [] ]
            }
            .set { ch_manta_input }

        MANTA_SOMATIC (
            ch_manta_input,
            ch_fasta,
            ch_fasta_fai,
            []  // config
        )
    }

    // =============================================
    // SOMATIC SNV/INDEL CALLING
    // =============================================

    if (doWF_SNV) {

        //
        // MODULE: Mutect2
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, [tbam, nbam], [tbai, nbai], [] ]
            }
            .set { ch_mutect2_input }

        // Prepare fai+gzi tuple
        ch_fasta_fai
            .map { meta, fai -> [ meta, fai, [] ] }
            .set { ch_fai_gzi }

        GATK4_MUTECT2 (
            ch_mutect2_input,
            ch_fasta,
            ch_fai_gzi,
            ch_dict,
            [],  // alleles
            [],  // alleles_tbi
            ch_germline_resource.map{ it[1] },
            ch_germline_resource_tbi.map{ it[1] },
            ch_pon.map{ it[1] },
            ch_pon_tbi.map{ it[1] }
        )

        //
        // MODULE: GetPileupSummaries (tumor)
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, tbam, tbai, [] ]
            }
            .set { ch_pileup_tumor_input }

        GETPILEUPSUMMARIES_TUMOR (
            ch_pileup_tumor_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_germline_resource.map{ it[1] },
            ch_germline_resource_tbi.map{ it[1] }
        )

        //
        // MODULE: GetPileupSummaries (normal)
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, nbam, nbai, [] ]
            }
            .set { ch_pileup_normal_input }

        GETPILEUPSUMMARIES_NORMAL (
            ch_pileup_normal_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_germline_resource.map{ it[1] },
            ch_germline_resource_tbi.map{ it[1] }
        )

        //
        // MODULE: CalculateContamination
        //
        GETPILEUPSUMMARIES_TUMOR.out.table
            .join(GETPILEUPSUMMARIES_NORMAL.out.table)
            .map { meta, tumor_table, normal_table ->
                [ meta, tumor_table, normal_table ]
            }
            .set { ch_contamination_input }

        GATK4_CALCULATECONTAMINATION ( ch_contamination_input )

        //
        // MODULE: LearnReadOrientationModel
        //
        GATK4_LEARNREADORIENTATIONMODEL (
            GATK4_MUTECT2.out.f1r2
        )

        //
        // MODULE: FilterMutectCalls
        //
        GATK4_MUTECT2.out.vcf
            .join(GATK4_MUTECT2.out.tbi)
            .join(GATK4_MUTECT2.out.stats)
            .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)
            .join(GATK4_CALCULATECONTAMINATION.out.contamination)
            .join(GATK4_CALCULATECONTAMINATION.out.segmentation)
            .map { meta, vcf, tbi, stats, orientationbias, contamination, segmentation ->
                [ meta, vcf, tbi, stats, orientationbias, segmentation, contamination, [] ]
            }
            .set { ch_filtermutect_input }

        GATK4_FILTERMUTECTCALLS (
            ch_filtermutect_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict
        )

        //
        // MODULE: Strelka2 somatic (uses Manta indel candidates)
        //
        ch_tumor_normal_pair
            .join(MANTA_SOMATIC.out.candidate_small_indels_vcf)
            .join(MANTA_SOMATIC.out.candidate_small_indels_vcf_tbi)
            .map { meta, tbam, tbai, nbam, nbai, manta_vcf, manta_tbi ->
                [ meta, nbam, nbai, tbam, tbai, manta_vcf, manta_tbi, [], [] ]
            }
            .set { ch_strelka_input }

        STRELKA_SOMATIC (
            ch_strelka_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] }
        )

        //
        // MODULE: Combine Strelka2 SNVs + Indels into single VCF
        //
        STRELKA_SOMATIC.out.vcf_snvs
            .join(STRELKA_SOMATIC.out.vcf_snvs_tbi)
            .join(STRELKA_SOMATIC.out.vcf_indels)
            .join(STRELKA_SOMATIC.out.vcf_indels_tbi)
            .map { meta, snv_vcf, snv_tbi, indel_vcf, indel_tbi ->
                [ meta, snv_vcf, snv_tbi, indel_vcf, indel_tbi ]
            }
            .set { ch_strelka_combine_input }

        STRELKA2_COMBINE_SOMATIC (
            ch_strelka_combine_input,
            ch_fasta,
            ch_fasta_fai
        )

        //
        // MODULE: Somatic Combine Channel (Mutect2 + Strelka2 union merge with annotations)
        //
        GATK4_FILTERMUTECTCALLS.out.vcf
            .join(GATK4_FILTERMUTECTCALLS.out.tbi)
            .join(
                ch_tumor_normal_pair.map { meta, tbam, tbai, nbam, nbai ->
                    [ meta, tbam, tbai, nbam, nbai ]
                }
            )
            .join(STRELKA2_COMBINE_SOMATIC.out.vcf.map { meta, vcf, tbi -> [ meta, vcf, tbi ] })
            .map { meta, mutect_vcf, mutect_tbi, tbam, tbai, nbam, nbai, strelka_vcf, strelka_tbi ->
                [ meta, mutect_vcf, mutect_tbi, tbam, tbai, nbam, nbai, strelka_vcf, strelka_tbi ]
            }
            .set { ch_somatic_combine_input }

        // Reference annotation files (optional, from params)
        ch_repeat_masker = params.repeat_masker
            ? Channel.value(file(params.repeat_masker, checkIfExists: true))
            : Channel.value([])
        ch_repeat_masker_tbi = params.repeat_masker_tbi
            ? Channel.value(file(params.repeat_masker_tbi, checkIfExists: true))
            : Channel.value([])
        ch_mapability_blacklist = params.mapability_blacklist
            ? Channel.value(file(params.mapability_blacklist, checkIfExists: true))
            : Channel.value([])
        ch_mapability_blacklist_tbi = params.mapability_blacklist_tbi
            ? Channel.value(file(params.mapability_blacklist_tbi, checkIfExists: true))
            : Channel.value([])
        ch_somatic_pon = params.somatic_pon
            ? Channel.value(file(params.somatic_pon, checkIfExists: true))
            : Channel.value([])
        ch_somatic_pon_tbi = params.somatic_pon_tbi
            ? Channel.value(file(params.somatic_pon_tbi, checkIfExists: true))
            : Channel.value([])
        ch_gnomad_somatic = params.gnomad_somatic
            ? Channel.value(file(params.gnomad_somatic, checkIfExists: true))
            : Channel.value([])
        ch_gnomad_somatic_tbi = params.gnomad_somatic_tbi
            ? Channel.value(file(params.gnomad_somatic_tbi, checkIfExists: true))
            : Channel.value([])

        SOMATIC_COMBINE_CHANNEL (
            ch_somatic_combine_input,
            ch_fasta,
            ch_fasta_fai,
            ch_repeat_masker,
            ch_repeat_masker_tbi,
            ch_mapability_blacklist,
            ch_mapability_blacklist_tbi,
            ch_somatic_pon,
            ch_somatic_pon_tbi,
            ch_gnomad_somatic,
            ch_gnomad_somatic_tbi
        )
        ch_versions = ch_versions.mix(SOMATIC_COMBINE_CHANNEL.out.versions.first())

        //
        // MODULE: Ensembl VEP annotation
        //
        GATK4_FILTERMUTECTCALLS.out.vcf
            .map { meta, vcf -> [ meta, vcf, [] ] }
            .set { ch_vep_input }

        ENSEMBLVEP_VEP (
            ch_vep_input,
            params.genome       ?: 'GRCh37',
            params.species      ?: 'homo_sapiens',
            params.vep_cache_version ?: '110',
            params.vep_cache    ? Channel.value(file(params.vep_cache, checkIfExists: true)) : Channel.value([]),
            ch_fasta,
            []  // extra_files
        )

        //
        // MODULE: vcf2maf
        //
        VCF2MAF (
            ENSEMBLVEP_VEP.out.vcf,
            ch_fasta.map{ it[1] },
            params.vep_cache ? Channel.value(file(params.vep_cache, checkIfExists: true)) : Channel.value([])
        )
        ch_versions = ch_versions.mix(VCF2MAF.out.versions.first())
    }

    // =============================================
    // SOMATIC SV CALLING
    // =============================================

    if (doWF_SV) {
        //
        // MODULE: Delly somatic SV calling (split by SV type, then merge)
        // Calls each SV type separately for parallelism and fault tolerance,
        // applies somatic + read support filtering, then merges per sample pair
        //

        // SV types to call in parallel
        ch_sv_types = Channel.from("DEL", "DUP", "INV", "BND", "INS")

        // Delly exclude regions (optional)
        ch_delly_exclude = params.delly_exclude_regions
            ? Channel.value(file(params.delly_exclude_regions, checkIfExists: true))
            : Channel.value([])

        DELLY_CALL_SOMATIC (
            ch_tumor_normal_pair,
            ch_sv_types,
            ch_fasta,
            ch_fasta_fai,
            ch_delly_exclude
        )
        ch_versions = ch_versions.mix(DELLY_CALL_SOMATIC.out.versions.first())

        //
        // MODULE: Delly combine - merge all SV types per sample pair
        //
        DELLY_CALL_SOMATIC.out.vcf
            .groupTuple(by: 0)
            .map { meta, svTypes, vcfs, tbis ->
                [ meta, vcfs, tbis ]
            }
            .set { ch_delly_combine_input }

        DELLY_COMBINE ( ch_delly_combine_input )
        ch_versions = ch_versions.mix(DELLY_COMBINE.out.versions.first())

        //
        // MODULE: Merge Delly + Manta somatic SVs
        // Concatenates Delly and Manta VCFs, filters PASS on canonical chromosomes
        //
        if (doWF_manta) {
            // Merge Delly + Manta somatic SVs (requires Manta output)
            DELLY_COMBINE.out.vcf
                .join(MANTA_SOMATIC.out.diploid_sv_vcf)
                .join(MANTA_SOMATIC.out.diploid_sv_vcf_tbi)
                .map { meta, delly_vcf, delly_tbi, manta_vcf, manta_tbi ->
                    [ meta, delly_vcf, delly_tbi, manta_vcf, manta_tbi ]
                }
                .set { ch_sv_merge_input }

            SOMATIC_MERGE_SV ( ch_sv_merge_input )
            ch_versions = ch_versions.mix(SOMATIC_MERGE_SV.out.versions.first())
        }

        //
        // MODULE: SvABA somatic SV calling
        //
        SVABA_SOMATIC (
            ch_tumor_normal_pair,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        //
        // MODULE: VCF2BEDPE somatic (Manta+Delly merged VCF → BEDPE)
        //
        if (doWF_manta) {
            SVTOOLS_VCF2BEDPE_SOMATIC (
                SOMATIC_MERGE_SV.out.vcf
            )

            //
            // MODULE: iAnnotateSV somatic (annotate BEDPE with blacklists + gene annotations)
            //
            ch_splice_sites = params.splice_sites
                ? Channel.value(file(params.splice_sites, checkIfExists: true))
                : Channel.value([])
            ch_sv_blacklist_bed = params.sv_blacklist_bed
                ? Channel.value(file(params.sv_blacklist_bed, checkIfExists: true))
                : Channel.value([])
            ch_sv_blacklist_bedpe = params.sv_blacklist_bedpe
                ? Channel.value(file(params.sv_blacklist_bedpe, checkIfExists: true))
                : Channel.value([])
            ch_sv_blacklist_foldback_bedpe = params.sv_blacklist_foldback_bedpe
                ? Channel.value(file(params.sv_blacklist_foldback_bedpe, checkIfExists: true))
                : Channel.value([])
            ch_sv_blacklist_te_bedpe = params.sv_blacklist_te_bedpe
                ? Channel.value(file(params.sv_blacklist_te_bedpe, checkIfExists: true))
                : Channel.value([])

            IANNOTATESV_SOMATIC (
                SVTOOLS_VCF2BEDPE_SOMATIC.out.bedpe,
                ch_repeat_masker,
                ch_mapability_blacklist,
                ch_sv_blacklist_bed,
                ch_sv_blacklist_bedpe,
                ch_sv_blacklist_foldback_bedpe,
                ch_sv_blacklist_te_bedpe,
                ch_splice_sites,
                params.genome ?: 'GRCh37'
            )

            //
            // MODULE: ClusterSV (cluster breakpoints)
            //
            CLUSTERSV ( IANNOTATESV_SOMATIC.out.bedpe_pass, params.genome ?: 'GRCh37' )

        }
    }

    // =============================================
    // FACETS (Copy Number)
    // =============================================

    if (doWF_facets) {
        //
        // MODULE: SNP-Pileup (local)
        // input: tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
        //        path facets_vcf
        //
        SNPPILEUP (
            ch_tumor_normal_pair,
            params.facets_vcf ? Channel.value(file(params.facets_vcf, checkIfExists: true)) : Channel.value([])
        )
        ch_versions = ch_versions.mix(SNPPILEUP.out.versions.first())

        //
        // MODULE: FACETS (local)
        //
        FACETS ( SNPPILEUP.out.pileup )
        ch_versions = ch_versions.mix(FACETS.out.versions.first())
    }

    // =============================================
    // FACETS PREVIEW QC & ANNOTATION
    // =============================================

    if (doWF_facets) {
        //
        // MODULE: FACETS Preview QC
        //
        FACETS_PREVIEW_QC (
            FACETS.out.facets_output
        )

        //
        // MODULE: Somatic FACETS Annotation (annotate MAF with CN + zygosity)
        // Requires both FACETS hisens_rdata and somatic MAF from VCF2MAF
        //
        if (doWF_SNV) {
            FACETS.out.hisens_rdata
                .join(VCF2MAF.out.maf)
                .set { ch_somatic_facets_anno_input }

            SOMATIC_FACETS_ANNOTATION (
                ch_somatic_facets_anno_input
            )
        }

    }

    // =============================================
    // SVCircos (circos plot visualization)
    // Requires both SV annotated BEDPE and FACETS CNV output
    // =============================================

    if (doWF_SV && doWF_facets) {
        IANNOTATESV_SOMATIC.out.bedpe_pass
            .join(FACETS.out.hisens_seg)
            .set { ch_svcircos_input }

        SVCIRCOS ( ch_svcircos_input, params.genome ?: 'GRCh37' )
    }

    // =============================================
    // MUTATION SIGNATURES (moved before blocks that depend on LOH/MSI)
    // =============================================

    if (doWF_mutSig) {
        //
        // MODULE: Mutational Signatures (tempoSig)
        //
        if (doWF_facets && doWF_SNV) {
            MUTSIG ( SOMATIC_FACETS_ANNOTATION.out.final_maf )
        }
    }

    // =============================================
    // SVCLONE
    // =============================================

    if (doWF_SV && doWF_facets && doWF_SNV && doWF_manta) {
        // SVclone needs: tumor/normal BAMs, annotated BEDPE, somatic MAF, FACETS CNV + ploidy
        ch_tumor_normal_pair
            .join(IANNOTATESV_SOMATIC.out.bedpe_pass)
            .join(SOMATIC_FACETS_ANNOTATION.out.final_maf)
            .join(FACETS.out.hisens_seg)
            .join(FACETS.out.purity)
            .set { ch_svclone_input }

        SVCLONE ( ch_svclone_input )
    }

    // =============================================
    // MSI
    // =============================================

    if (doWF_msiSensor) {
        //
        // MODULE: MSIsensor-pro scan
        // input: tuple val(meta), path(fasta)
        //
        MSISENSORPRO_SCAN ( ch_fasta )

        //
        // MODULE: MSIsensor-pro msi somatic
        // input: tuple val(meta), path(normal), path(normal_index),
        //              path(tumor), path(tumor_index), path(intervals)
        //        tuple val(meta2), path(fasta)
        //        path(msisensor_scan)
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, nbam, nbai, tbam, tbai, [] ]
            }
            .set { ch_msi_input }

        MSISENSORPRO_MSISOMATIC (
            ch_msi_input,
            ch_fasta,
            MSISENSORPRO_SCAN.out.list.map{ it[1] }
        )
    }

    // =============================================
    // HLA TYPING & LOH
    // =============================================

    if (doWF_loh) {
        //
        // MODULE: Polysolver (local) - HLA typing on normal BAMs
        // input: tuple val(meta), path(bam), path(bai)
        //
        POLYSOLVER (
            ch_recal_branched.normal.map { meta, bam, bai -> [ meta, bam, bai ] }
        )
        ch_versions = ch_versions.mix(POLYSOLVER.out.versions.first())

        //
        // MODULE: LOHHLA (local) - HLA LOH detection
        // input: tuple val(meta), path(tumor_bam), path(tumor_bai),
        //              path(normal_bam), path(normal_bai), path(hla_types)
        //        path hla_fasta
        //        path hla_dat
        //
        // Pair tumor-normal with Polysolver HLA types from normal
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta.normal_id, meta, tbam, tbai, nbam, nbai ]
            }
            .combine(
                POLYSOLVER.out.hla_types.map { meta, hla -> [ meta.sample ?: meta.id, hla ] },
                by: 0
            )
            .map { normal_id, meta, tbam, tbai, nbam, nbai, hla_types ->
                [ meta, tbam, tbai, nbam, nbai, hla_types ]
            }
            .set { ch_lohhla_input }

        LOHHLA (
            ch_lohhla_input,
            params.hla_fasta ? Channel.value(file(params.hla_fasta, checkIfExists: true)) : Channel.value([]),
            params.hla_dat   ? Channel.value(file(params.hla_dat, checkIfExists: true))   : Channel.value([])
        )
        ch_versions = ch_versions.mix(LOHHLA.out.versions.first())
    }

    // =============================================
    // GERMLINE VARIANT CALLING
    // =============================================

    if (doWF_germSNV) {
        //
        // MODULE: GATK4 HaplotypeCaller
        // input: tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
        //        tuple val(meta2), path(fasta)
        //        tuple val(meta3), path(fai)
        //        tuple val(meta4), path(dict)
        //        tuple val(meta5), path(dbsnp)
        //        tuple val(meta6), path(dbsnp_tbi)
        //
        ch_recal_branched.normal
            .map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }
            .set { ch_hc_input }

        GATK4_HAPLOTYPECALLER (
            ch_hc_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi
        )

        //
        // MODULE: Strelka2 germline
        // input: tuple val(meta), path(bam), path(bai)
        //        tuple val(meta2), path(fasta)
        //        tuple val(meta3), path(fai)
        //        path(call_regions), path(call_regions_tbi)
        //
        ch_call_regions = params.call_regions
            ? Channel.value(file(params.call_regions, checkIfExists: true))
            : Channel.value([])
        ch_call_regions_tbi = params.call_regions_tbi
            ? Channel.value(file(params.call_regions_tbi, checkIfExists: true))
            : Channel.value([])

        STRELKA2_GERMLINE (
            ch_recal_branched.normal,
            ch_fasta,
            ch_fasta_fai,
            ch_call_regions,
            ch_call_regions_tbi
        )
        ch_versions = ch_versions.mix(STRELKA2_GERMLINE.out.versions.first())

        //
        // MODULE: Germline Combine Channel (HaplotypeCaller + Strelka2 union merge)
        //
        //
        // MODULE: Germline Combine Channel (HaplotypeCaller + Strelka2 union merge)
        //
        // Join HC output with Strelka2 germline output
        GATK4_HAPLOTYPECALLER.out.vcf
            .join(GATK4_HAPLOTYPECALLER.out.tbi)
            .join(STRELKA2_GERMLINE.out.vcf.map { meta, vcf, tbi -> [ meta, vcf, tbi ] })
            .map { meta, hc_vcf, hc_tbi, strelka_vcf, strelka_tbi ->
                [ meta.patient, meta, hc_vcf, hc_tbi, strelka_vcf, strelka_tbi ]
            }
            .combine(
                ch_recal_branched.tumor.map { meta, bam, bai -> [ meta.patient, meta, bam, bai ] },
                by: 0
            )
            .map { patient, normal_meta, hc_vcf, hc_tbi, strelka_vcf, strelka_tbi, tumor_meta, tbam, tbai ->
                def pair_meta = [
                    id:        "${tumor_meta.sample}__${normal_meta.sample}",
                    patient:   patient,
                    tumor_id:  tumor_meta.sample,
                    normal_id: normal_meta.sample
                ]
                [ pair_meta, hc_vcf, hc_tbi, strelka_vcf, strelka_tbi, tbam, tbai ]
            }
            .set { ch_germline_combine_input }

        // Reference annotation files for germline
        ch_repeat_masker_germ = params.repeat_masker
            ? Channel.value(file(params.repeat_masker, checkIfExists: true))
            : Channel.value([])
        ch_repeat_masker_tbi_germ = params.repeat_masker_tbi
            ? Channel.value(file(params.repeat_masker_tbi, checkIfExists: true))
            : Channel.value([])
        ch_mapability_blacklist_germ = params.mapability_blacklist
            ? Channel.value(file(params.mapability_blacklist, checkIfExists: true))
            : Channel.value([])
        ch_mapability_blacklist_tbi_germ = params.mapability_blacklist_tbi
            ? Channel.value(file(params.mapability_blacklist_tbi, checkIfExists: true))
            : Channel.value([])
        ch_gnomad_germline = params.gnomad_germline
            ? Channel.value(file(params.gnomad_germline, checkIfExists: true))
            : Channel.value([])
        ch_gnomad_germline_tbi = params.gnomad_germline_tbi
            ? Channel.value(file(params.gnomad_germline_tbi, checkIfExists: true))
            : Channel.value([])

        GERMLINE_COMBINE_CHANNEL (
            ch_germline_combine_input,
            ch_fasta,
            ch_fasta_fai,
            ch_repeat_masker_germ,
            ch_repeat_masker_tbi_germ,
            ch_mapability_blacklist_germ,
            ch_mapability_blacklist_tbi_germ,
            ch_gnomad_germline,
            ch_gnomad_germline_tbi
        )
        ch_versions = ch_versions.mix(GERMLINE_COMBINE_CHANNEL.out.versions.first())

        //
        // MODULE: Germline VEP + MAF annotation
        //
        ch_germline_vep_cache = params.vep_cache
            ? Channel.value(file(params.vep_cache, checkIfExists: true))
            : Channel.value([])
        ch_germline_isoforms = params.vep_custom_isoforms
            ? Channel.value(file(params.vep_custom_isoforms, checkIfExists: true))
            : Channel.value([])

        GERMLINE_ANNOTATE_MAF (
            GERMLINE_COMBINE_CHANNEL.out.germline_vcf,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_germline_vep_cache,
            ch_germline_isoforms
        )
        ch_versions = ch_versions.mix(GERMLINE_ANNOTATE_MAF.out.versions.first())
    }

    // =============================================
    // GERMLINE SV CALLING
    // =============================================

    if (doWF_germSV) {
        //
        // MODULE: Delly germline SV calling (split by SV type)
        //
        ch_sv_types_germline = Channel.from("DEL", "DUP", "INV", "BND", "INS")

        ch_delly_exclude_germline = params.delly_exclude_regions
            ? Channel.value(file(params.delly_exclude_regions, checkIfExists: true))
            : Channel.value([])

        DELLY_CALL_GERMLINE (
            ch_recal_branched.normal,
            ch_sv_types_germline,
            ch_fasta,
            ch_fasta_fai,
            ch_delly_exclude_germline
        )
        ch_versions = ch_versions.mix(DELLY_CALL_GERMLINE.out.versions.first())

        //
        // MODULE: Manta germline
        //
        ch_manta_germline_regions = params.sv_calling_include_regions
            ? Channel.value(file(params.sv_calling_include_regions, checkIfExists: true))
            : Channel.value([])
        ch_manta_germline_regions_tbi = params.sv_calling_include_regions_tbi
            ? Channel.value(file(params.sv_calling_include_regions_tbi, checkIfExists: true))
            : Channel.value([])

        MANTA_GERMLINE (
            ch_recal_branched.normal,
            ch_fasta,
            ch_fasta_fai,
            ch_manta_germline_regions,
            ch_manta_germline_regions_tbi
        )
        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())

        //
        // MODULE: Merge Delly + Manta germline SVs
        //
        DELLY_CALL_GERMLINE.out.vcf
            .groupTuple(by: 0)
            .map { meta, svTypes, vcfs, tbis -> [ meta, vcfs, tbis ] }
            .combine(MANTA_GERMLINE.out.vcf.map { meta, vcf, tbi -> [ meta, vcf, tbi ] }, by: 0)
            .map { meta, delly_vcfs, delly_tbis, manta_vcf, manta_tbi ->
                [ meta, delly_vcfs, delly_tbis, manta_vcf, manta_tbi ]
            }
            .set { ch_germline_sv_merge_input }

        GERMLINE_MERGE_SV ( ch_germline_sv_merge_input )
        ch_versions = ch_versions.mix(GERMLINE_MERGE_SV.out.versions.first())

        //
        // MODULE: SvABA germline SV calling
        //
        SVABA_GERMLINE (
            ch_recal_branched.normal,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        //
        // MODULE: VCF2BEDPE germline
        //
        SVTOOLS_VCF2BEDPE_GERMLINE (
            GERMLINE_MERGE_SV.out.vcf
        )

        //
        // MODULE: iAnnotateSV germline
        //
        ch_splice_sites_germ = params.splice_sites
            ? Channel.value(file(params.splice_sites, checkIfExists: true))
            : Channel.value([])
        ch_sv_blacklist_bed_germ = params.sv_blacklist_bed
            ? Channel.value(file(params.sv_blacklist_bed, checkIfExists: true))
            : Channel.value([])
        ch_sv_blacklist_bedpe_germ = params.sv_blacklist_bedpe
            ? Channel.value(file(params.sv_blacklist_bedpe, checkIfExists: true))
            : Channel.value([])
        ch_sv_blacklist_foldback_bedpe_germ = params.sv_blacklist_foldback_bedpe
            ? Channel.value(file(params.sv_blacklist_foldback_bedpe, checkIfExists: true))
            : Channel.value([])
        ch_sv_blacklist_te_bedpe_germ = params.sv_blacklist_te_bedpe
            ? Channel.value(file(params.sv_blacklist_te_bedpe, checkIfExists: true))
            : Channel.value([])
        ch_repeat_masker_sv_germ = params.repeat_masker
            ? Channel.value(file(params.repeat_masker, checkIfExists: true))
            : Channel.value([])
        ch_mapability_blacklist_sv_germ = params.mapability_blacklist
            ? Channel.value(file(params.mapability_blacklist, checkIfExists: true))
            : Channel.value([])

        IANNOTATESV_GERMLINE (
            SVTOOLS_VCF2BEDPE_GERMLINE.out.bedpe,
            ch_repeat_masker_sv_germ,
            ch_mapability_blacklist_sv_germ,
            ch_sv_blacklist_bed_germ,
            ch_sv_blacklist_bedpe_germ,
            ch_sv_blacklist_foldback_bedpe_germ,
            ch_sv_blacklist_te_bedpe_germ,
            ch_splice_sites_germ,
            params.genome ?: 'GRCh37'
        )
    }

    // =============================================
    // GERMLINE FACETS ANNOTATION (after germline SNV so GERMLINE_COMBINE_CHANNEL is available)
    // =============================================

    if (doWF_facets && doWF_germSNV) {
        FACETS.out.hisens_rdata
            .join(GERMLINE_COMBINE_CHANNEL.out.germline_vcf)
            .set { ch_germline_facets_anno_input }

        GERMLINE_FACETS_ANNOTATION (
            ch_germline_facets_anno_input
        )
    }

    // =============================================
    // NEOANTIGEN PREDICTION (after LOH/POLYSOLVER so POLYSOLVER.out is available)
    // =============================================

    if (doWF_loh && doWF_SNV) {
        ch_neoantigen_cdna = params.neoantigen_cdna
            ? Channel.value(file(params.neoantigen_cdna, checkIfExists: true))
            : Channel.value([])
        ch_neoantigen_cds = params.neoantigen_cds
            ? Channel.value(file(params.neoantigen_cds, checkIfExists: true))
            : Channel.value([])

        // Pair polysolver output with somatic MAF
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai -> [ meta.normal_id, meta ] }
            .combine(
                POLYSOLVER.out.hla_types.map { meta, hla -> [ meta.sample ?: meta.id, hla ] },
                by: 0
            )
            .map { normal_id, meta, hla -> [ meta.id, meta, hla ] }
            .combine(
                VCF2MAF.out.maf.map { meta, maf -> [ meta.id, maf ] },
                by: 0
            )
            .map { id, meta, hla, maf -> [ meta, hla, maf ] }
            .set { ch_neoantigen_input }

        NEOANTIGEN (
            ch_neoantigen_input,
            ch_neoantigen_cdna,
            ch_neoantigen_cds
        )
    }

    // =============================================
    // METADATA PARSER (after MSI + LOH + MUTSIG so all inputs available)
    // =============================================

    if (doWF_mdParse) {
        FACETS.out.purity
            .join(SOMATIC_FACETS_ANNOTATION.out.final_maf)
            .join(FACETS_PREVIEW_QC.out.facets_qc)
            .join(MSISENSORPRO_MSISOMATIC.out.output_report)
            .join(MUTSIG.out.mutsig_results)
            .set { ch_metadata_partial }

        ch_metadata_partial
            .map { meta, purity, maf, qc, msi, mutsig ->
                [ meta.normal_id, meta, purity, maf, qc, msi, mutsig ]
            }
            .combine(
                POLYSOLVER.out.hla_types.map { meta, hla -> [ meta.sample ?: meta.id, hla ] },
                by: 0
            )
            .map { nid, meta, purity, maf, qc, msi, mutsig, hla ->
                [ meta, purity, maf, qc, msi, mutsig, hla ]
            }
            .set { ch_metadata_input }

        ch_coding_bed = params.coding_bed
            ? Channel.value(file(params.coding_bed, checkIfExists: true))
            : Channel.value([])

        METADATA_PARSER (
            ch_metadata_input,
            ch_coding_bed
        )
    }

    // =============================================
    // WGS-ONLY: ASCAT + BRASS + HRDetect
    // =============================================

    if (isWGS && doWF_SV) {
        //
        // MODULE: ASCAT AlleleCount (Somatic CNV calling - WGS only)
        //
        ASCAT_ALLELECOUNT (
            ch_tumor_normal_pair,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.snp_gc_corrections ? Channel.value(file(params.snp_gc_corrections, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        //
        // MODULE: ASCAT Run (Somatic CNV calling - WGS only)
        //
        ASCAT_ALLELECOUNT.out.alleles
            .groupTuple(by: 0)
            .map { meta, allele_counts ->
                def ascat_tar = allele_counts.size() > 0 ? allele_counts[0] : file('NO_FILE')
                [ meta, ascat_tar ]
            }
            .join(ch_tumor_normal_pair)
            .map { meta, ascat_tar, tbam, tbai, nbam, nbai ->
                [ meta, ascat_tar, tbam, tbai, nbam, nbai ]
            }
            .set { ch_ascat_input }

        ASCAT_RUN (
            ch_ascat_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.snp_gc_corrections ? Channel.value(file(params.snp_gc_corrections, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        //
        // MODULE: BRASS GenerateBas (Somatic rearrangement annotation - WGS only)
        // Run separately on tumor and normal BAMs
        //
        ch_recal_bam_bai
            .map { meta, bam, bai ->
                def bas_meta = [
                    id:     "${meta.sample}__bas",
                    sample: meta.sample,
                    status: meta.status
                ]
                [ bas_meta, bam, bai ]
            }
            .set { ch_brass_bas_input }

        BRASS_GENERATE_BAS (
            ch_brass_bas_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] }
        )

        //
        // MODULE: BRASS Input (Somatic rearrangement annotation - WGS only)
        // Requires tumor + normal BAS files
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai -> [ meta.tumor_id, meta, tbam, tbai, nbam, nbai ] }
            .combine(
                BRASS_GENERATE_BAS.out.bas.map { meta, bas ->
                    [ meta.sample, meta, bas ]
                },
                by: 0
            )
            .map { tid, meta, tbam, tbai, nbam, nbai, tumor_meta, tumor_bas ->
                [ meta.normal_id, meta, tbam, tbai, nbam, nbai, tumor_bas ]
            }
            .combine(
                BRASS_GENERATE_BAS.out.bas.map { meta, bas ->
                    [ meta.sample, bas ]
                },
                by: 0
            )
            .map { nid, meta, tbam, tbai, nbam, nbai, tumor_bas, normal_bas ->
                [ meta, tbam, tbai, tumor_bas, nbam, nbai, normal_bas ]
            }
            .set { ch_brass_input_input }

        BRASS_INPUT (
            ch_brass_input_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.brass_ref_dir ? Channel.value(file(params.brass_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE')),
            params.vagrent_ref_dir ? Channel.value(file(params.vagrent_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        //
        // MODULE: BRASS Cover (Somatic rearrangement annotation - WGS only)
        //
        BRASS_COVER (
            ch_brass_input_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.brass_ref_dir ? Channel.value(file(params.brass_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE')),
            params.vagrent_ref_dir ? Channel.value(file(params.vagrent_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        //
        // MODULE: BRASS Run (Somatic rearrangement annotation - WGS only)
        // Combines BRASS input + cover outputs with BAMs + ASCAT results
        //
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai -> [ meta.tumor_id, meta, tbam, tbai, nbam, nbai ] }
            .combine(
                BRASS_GENERATE_BAS.out.bas.map { meta, bas ->
                    [ meta.sample, meta, bas ]
                },
                by: 0
            )
            .map { tid, meta, tbam, tbai, nbam, nbai, tumor_meta, tumor_bas ->
                [ meta.normal_id, meta, tbam, tbai, nbam, nbai, tumor_bas ]
            }
            .combine(
                BRASS_GENERATE_BAS.out.bas.map { meta, bas ->
                    [ meta.sample, bas ]
                },
                by: 0
            )
            .map { nid, meta, tbam, tbai, nbam, nbai, tumor_bas, normal_bas ->
                [ meta, tbam, tbai, tumor_bas, nbam, nbai, normal_bas ]
            }
            .join(BRASS_INPUT.out.brass_input)
            .join(BRASS_COVER.out.brass_cover)
            .join(ASCAT_RUN.out.ascat_samplestatistics)
            .map { meta, tbam, tbai, tumor_bas, nbam, nbai, normal_bas, brass_input_data, brass_cover_data, ascat_stats ->
                [ meta, tbam, tbai, tumor_bas, nbam, nbai, normal_bas, brass_input_data, brass_cover_data, ascat_stats ]
            }
            .set { ch_brass_run_input }

        BRASS_RUN (
            ch_brass_run_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.brass_ref_dir ? Channel.value(file(params.brass_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE')),
            params.vagrent_ref_dir ? Channel.value(file(params.vagrent_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        //
        // MODULE: HRDetect (Homologous recombination deficiency detection - WGS only)
        // Requires: somatic MAF + FACETS CNV + annotated BEDPE + HRdetect script
        //
        if (doWF_SNV && doWF_facets && doWF_manta) {
            ch_hrdetect_script = params.hrdetect_script
                ? Channel.value(file(params.hrdetect_script, checkIfExists: true))
                : Channel.value(file('NO_FILE'))

            SOMATIC_FACETS_ANNOTATION.out.final_maf
                .join(FACETS.out.hisens_seg)
                .join(IANNOTATESV_SOMATIC.out.bedpe_pass)
                .map { meta, maf, cnv, bedpe ->
                    [ meta, maf, cnv, bedpe ]
                }
                .set { ch_hrdetect_input }

            HRDETECT (
                ch_hrdetect_input,
                ch_hrdetect_script
            )
        }

        //
        // MODULE: SV Signatures (Somatic rearrangement signatures - WGS only)
        // Requires: annotated BEDPE + SV signature script
        //
        if (doWF_manta) {
            ch_svsig_script = params.sv_signature_script
                ? Channel.value(file(params.sv_signature_script, checkIfExists: true))
                : Channel.value(file('NO_FILE'))

            IANNOTATESV_SOMATIC.out.bedpe_pass
                .set { ch_svsig_input }

            SV_SIGNATURES (
                ch_svsig_input,
                ch_svsig_script
            )
        }
    }

    // =============================================
    // QC
    // =============================================

    if (doWF_QC) {
        //
        // MODULE: Picard CollectHsMetrics (nf-core)
        // input: tuple val(meta), path(bam), path(bai), path(bait_intervals), path(target_intervals)
        //        tuple val(meta2), path(ref)
        //        tuple val(meta3), path(ref_fai)
        //        tuple val(meta4), path(ref_dict)
        //        tuple val(meta5), path(ref_gzi)
        //
        ch_recal_bam_bai
            .map { meta, bam, bai ->
                [ meta, bam, bai, params.bait_intervals ? file(params.bait_intervals) : [], params.target_intervals ? file(params.target_intervals) : [] ]
            }
            .set { ch_hsmetrics_input }

        PICARD_COLLECTHSMETRICS (
            ch_hsmetrics_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            Channel.value([ [:], [] ])   // ref_gzi (empty)
        )
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTHSMETRICS.out.metrics.collect{it[1]})

        //
        // MODULE: Qualimap BAM QC (nf-core)
        // input: tuple val(meta), path(bam)
        //        path gff
        //
        QUALIMAP_BAMQC (
            ch_recal_bam_bai.map { meta, bam, bai -> [ meta, bam ] },
            []  // gff
        )

        //
        // MODULE: Conpair pileup (local) - run on all recalibrated BAMs
        // input: tuple val(meta), path(bam), path(bai)
        //        path fasta
        //        path fasta_fai
        //        path dict
        //
        CONPAIR_PILEUP (
            ch_recal_bam_bai.map { meta, bam, bai -> [ meta, bam, bai ] },
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )
        ch_versions = ch_versions.mix(CONPAIR_PILEUP.out.versions.first())

        //
        // MODULE: Conpair concordance (local) - pair tumor/normal pileups
        // input: tuple val(meta), path(tumor_pileup), path(normal_pileup)
        //
        CONPAIR_PILEUP.out.pileup
            .branch {
                tumor:  it[0].status == 1
                normal: it[0].status == 0
            }
            .set { ch_conpair_branched }

        ch_conpair_branched.tumor
            .map { meta, pileup -> [ meta.patient, meta, pileup ] }
            .combine(
                ch_conpair_branched.normal.map { meta, pileup -> [ meta.patient, meta, pileup ] },
                by: 0
            )
            .map { patient, tumor_meta, tumor_pileup, normal_meta, normal_pileup ->
                def pair_meta = [
                    id:        "${tumor_meta.sample}__${normal_meta.sample}",
                    patient:   patient,
                    tumor_id:  tumor_meta.sample,
                    normal_id: normal_meta.sample
                ]
                [ pair_meta, tumor_pileup, normal_pileup ]
            }
            .set { ch_conpair_concordance_input }

        CONPAIR_CONCORDANCE ( ch_conpair_concordance_input )
        ch_versions = ch_versions.mix(CONPAIR_CONCORDANCE.out.versions.first())

        //
        // MODULE: Conpair All (combined concordance + contamination)
        //
        CONPAIR_ALL (
            ch_conpair_concordance_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        //
        // MODULE: Alfred QC (BAM quality metrics)
        //
        ch_recal_bam_bai
            .map { meta, bam, bai ->
                def targets = params.target_intervals ? file(params.target_intervals) : file('NO_FILE')
                def targets_idx = params.target_intervals ? file("${params.target_intervals}.idx", checkIfExists: false) : file('NO_FILE2')
                [ meta, bam, bai, targets, targets_idx ]
            }
            .set { ch_alfred_input }

        ALFRED (
            ch_alfred_input,
            ch_fasta.map{ it[1] }
        )

        //
        // MODULE: MultiQC Sample-level report
        // Per-sample QC metrics aggregation
        //
        ch_multiqc_sample_configs = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

        ALFRED.out.alfred_qc
            .map { meta, rg_n, rg_y ->
                [ meta.sample, meta, rg_n, rg_y ]
            }
            .combine(
                FASTP.out.json.map { meta, json -> [ meta.sample, json ] },
                by: 0
            )
            .map { sample, meta, rg_n, rg_y, fastp_json ->
                [ meta, rg_n, rg_y, fastp_json ]
            }
            .join(QUALIMAP_BAMQC.out.results)
            .map { meta, rg_n, rg_y, fastp_json, qualimap_dir ->
                [ meta, rg_n, rg_y, fastp_json, qualimap_dir ]
            }
            .join(PICARD_COLLECTHSMETRICS.out.metrics)
            .map { meta, rg_n, rg_y, fastp_json, qualimap_dir, hsmetrics ->
                [ meta, rg_n, rg_y, fastp_json, qualimap_dir, hsmetrics ]
            }
            .set { ch_multiqc_sample_input }

        MULTIQC_SAMPLE (
            ch_multiqc_sample_input,
            ch_multiqc_sample_configs.toList()
        )

        //
        // MODULE: MultiQC Somatic pair-level report
        // Per tumor-normal pair QC metrics aggregation
        //
        if (doWF_SNV && doWF_facets) {
            CONPAIR_CONCORDANCE.out.concordance
                .join(FACETS.out.summary)
                .join(FACETS_PREVIEW_QC.out.facets_qc)
                .map { meta, conpair, facets_sum, facets_qc ->
                    [ meta, conpair, facets_sum, facets_qc ]
                }
                .set { ch_multiqc_somatic_input }

            MULTIQC_SOMATIC (
                ch_multiqc_somatic_input,
                ch_multiqc_sample_configs.toList()
            )
        }
    }

    // =============================================
    // MultiQC (nf-core - fallback aggregation)
    // =============================================

    if (doWF_QC && !params.skip_multiqc) {

        ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)   : Channel.empty()

        //
        // MODULE: MultiQC (nf-core)
        // input: path multiqc_files, stageAs: "?/*"
        //        path(multiqc_config)
        //        path(extra_multiqc_config)
        //        path(multiqc_logo)
        //        path(replace_names)
        //        path(sample_names)
        //
        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],  // replace_names
            []   // sample_names
        )
    }

    // =============================================
    // AGGREGATION (cohort-level outputs)
    // =============================================

    if (params.aggregate) {
        //
        // MODULE: Aggregate Somatic MAF (union of all somatic MAF files)
        //
        if (doWF_SNV && doWF_facets) {
            SOMATIC_FACETS_ANNOTATION.out.final_maf
                .map { meta, maf -> maf }
                .collect()
                .set { ch_aggregate_somatic_maf_input }

            AGGREGATE_SOMATIC_MAF ( ch_aggregate_somatic_maf_input )
        }

        //
        // MODULE: Aggregate Somatic SV (union of all somatic SV BEDPE files)
        //
        if (doWF_SV && doWF_manta) {
            IANNOTATESV_SOMATIC.out.bedpe_pass
                .map { meta, bedpe -> bedpe }
                .collect()
                .set { ch_aggregate_somatic_sv_input }

            AGGREGATE_SOMATIC_SV ( ch_aggregate_somatic_sv_input )
        }

        //
        // MODULE: Aggregate Somatic FACETS (union of all FACETS CNV files)
        //
        if (doWF_facets) {
            FACETS.out.hisens_seg
                .map { meta, seg -> seg }
                .collect()
                .set { ch_aggregate_somatic_facets_input }

            AGGREGATE_SOMATIC_FACETS ( ch_aggregate_somatic_facets_input )
        }

        //
        // MODULE: Aggregate Somatic NetMHC (neoantigen predictions)
        //
        if (doWF_loh && doWF_SNV) {
            NEOANTIGEN.out.predictions
                .map { meta, neo -> neo }
                .collect()
                .set { ch_aggregate_somatic_netmhc_input }

            AGGREGATE_SOMATIC_NETMHC ( ch_aggregate_somatic_netmhc_input )
        }

        //
        // MODULE: Aggregate Somatic Metadata (clinical metadata summary)
        //
        if (doWF_mdParse) {
            METADATA_PARSER.out.metadata
                .map { meta, metadata -> metadata }
                .collect()
                .set { ch_aggregate_somatic_metadata_input }

            AGGREGATE_SOMATIC_METADATA ( ch_aggregate_somatic_metadata_input )
        }

        //
        // MODULE: Aggregate Somatic LOH/HLA (union of LOH HLA output)
        //
        if (doWF_loh) {
            LOHHLA.out.predictions
                .map { meta, summary -> summary }
                .collect()
                .set { ch_aggregate_somatic_lohhla_input }

            AGGREGATE_SOMATIC_LOHHLA ( ch_aggregate_somatic_lohhla_input )
        }

        //
        // MODULE: Aggregate Somatic HRDetect (homologous recombination deficiency)
        //
        if (isWGS && doWF_SV && doWF_SNV && doWF_facets) {
            HRDETECT.out.hrdetect_output
                .map { meta, hrdetect -> hrdetect }
                .collect()
                .set { ch_aggregate_somatic_hrdetect_input }

            AGGREGATE_SOMATIC_HRDETECT ( ch_aggregate_somatic_hrdetect_input )
        }

        //
        // MODULE: Aggregate Somatic SVClone (clonal SV analysis)
        //
        if (doWF_SV && doWF_facets && doWF_SNV && doWF_manta) {
            SVCLONE.out.cluster_certainty
                .map { meta, sv_cert, snv_cert -> [ sv_cert, snv_cert ] }
                .collect()
                .set { ch_aggregate_somatic_svclone_input }

            AGGREGATE_SOMATIC_SVCLONE ( ch_aggregate_somatic_svclone_input )
        }

        //
        // MODULE: Aggregate Somatic SV Signatures (rearrangement signatures)
        //
        if (isWGS && doWF_SV && doWF_manta) {
            SV_SIGNATURES.out.sv_signatures
                .map { meta, sigs -> sigs }
                .collect()
                .set { ch_aggregate_somatic_svsignatures_input }

            AGGREGATE_SOMATIC_SVSIGNATURES ( ch_aggregate_somatic_svsignatures_input )
        }

        //
        // MODULE: Aggregate Germline MAF (union of all germline MAF files)
        //
        if (doWF_germSNV) {
            GERMLINE_ANNOTATE_MAF.out.maf_file
                .map { meta, maf -> maf }
                .collect()
                .set { ch_aggregate_germline_maf_input }

            AGGREGATE_GERMLINE_MAF ( ch_aggregate_germline_maf_input )
        }

        //
        // MODULE: Aggregate Germline SV (union of all germline SV BEDPE files)
        //
        if (doWF_germSV) {
            IANNOTATESV_GERMLINE.out.bedpe_pass
                .map { meta, bedpe -> bedpe }
                .collect()
                .set { ch_aggregate_germline_sv_input }

            AGGREGATE_GERMLINE_SV ( ch_aggregate_germline_sv_input )
        }

        //
        // MODULE: Aggregate QC BAM (union of BAM quality metrics)
        //
        if (doWF_QC) {
            ALFRED.out.alfred_qc
                .map { meta, rg_n, rg_y -> [ rg_n, rg_y ] }
                .flatten()
                .collect()
                .set { ch_aggregate_alfred_input }

            PICARD_COLLECTHSMETRICS.out.metrics
                .map { meta, metrics -> metrics }
                .collect()
                .set { ch_aggregate_hsmetrics_input }

            AGGREGATE_QC_BAM ( ch_aggregate_alfred_input, ch_aggregate_hsmetrics_input )
        }

        //
        // MODULE: Aggregate QC Conpair (concordance metrics across cohort)
        //
        if (doWF_QC) {
            CONPAIR_CONCORDANCE.out.concordance
                .map { meta, concordance -> concordance }
                .collect()
                .set { ch_aggregate_qc_conpair_input }

            AGGREGATE_QC_CONPAIR ( ch_aggregate_qc_conpair_input )
        }
    }

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
