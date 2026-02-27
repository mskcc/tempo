/*
~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                       } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                      } from '../modules/nf-core/multiqc/main'
include { FASTP                                        } from '../modules/nf-core/fastp/main'
include { BWAMEM2_MEM                                  } from '../modules/nf-core/bwamem2/mem/main'
include { EXTRACT_READ_GROUP                             } from '../modules/local/extract_read_group/main'
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
include { GATK4_MERGEVCFS                              } from '../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_SOMATIC                              } from '../modules/nf-core/strelka/somatic/main'
include { MANTA_SOMATIC                                } from '../modules/nf-core/manta/somatic/main'
include { ENSEMBLVEP_VEP                               } from '../modules/nf-core/ensemblvep/vep/main'
include { MSISENSORPRO_SCAN                            } from '../modules/nf-core/msisensorpro/scan/main'
include { MSISENSORPRO_MSISOMATIC                      } from '../modules/nf-core/msisensorpro/msisomatic/main'
include { PICARD_COLLECTHSMETRICS                      } from '../modules/nf-core/picard/collecthsmetrics/main'
include { QUALIMAP_BAMQC                               } from '../modules/nf-core/qualimap/bamqc/main'

include { POLYSOLVER                  } from '../modules/local/polysolver/main'
include { LOHHLA                      } from '../modules/local/lohhla/main'
include { SNPPILEUP                   } from '../modules/local/snppileup/main'
include { FACETS                      } from '../modules/local/facets/main'
include { CONPAIR_PILEUP              } from '../modules/local/conpair/pileup/main'
include { CONPAIR_CONCORDANCE         } from '../modules/local/conpair/concordance/main'
include { VCF2MAF                     } from '../modules/local/vcf2maf/main'
include { DELLY_CALL_SOMATIC          } from '../modules/local/delly/call/main'
include { DELLY_COMBINE               } from '../modules/local/delly/combine/main'
include { DELLY_COMBINE as DELLY_COMBINE_GERMLINE } from '../modules/local/delly/combine/main'
include { DELLY_CALL_GERMLINE         } from '../modules/local/delly/call_germline/main'
include { STRELKA2_COMBINE_SOMATIC    } from '../modules/local/strelka2/combine_somatic/main'
include { SOMATIC_COMBINE_CHANNEL     } from '../modules/local/somatic/combine_channel/main'
include { SOMATIC_MERGE_SV            } from '../modules/local/somatic/merge_sv/main'
include { STRELKA2_GERMLINE           } from '../modules/local/strelka2/germline/main'
include { MANTA_GERMLINE              } from '../modules/local/manta/germline/main'
include { GERMLINE_MERGE_SV           } from '../modules/local/germline/merge_sv/main'
include { GERMLINE_COMBINE_CHANNEL    } from '../modules/local/germline/combine_channel/main'

// SV pipeline modules
include { SVABA_SOMATIC               } from '../modules/local/svaba/somatic/main'
include { SVABA_GERMLINE              } from '../modules/local/svaba/germline/main'
include { SVTOOLS_VCF2BEDPE_SOMATIC   } from '../modules/local/svtools/vcf2bedpe_somatic/main'
include { SVTOOLS_VCF2BEDPE_GERMLINE  } from '../modules/local/svtools/vcf2bedpe_germline/main'
include { IANNOTATESV_SOMATIC         } from '../modules/local/iannotatesv/somatic/main'
include { IANNOTATESV_GERMLINE        } from '../modules/local/iannotatesv/germline/main'
include { CLUSTERSV                   } from '../modules/local/clustersv/main'
include { SVCIRCOS                    } from '../modules/local/svcircos/main'
include { SVCLONE                     } from '../modules/local/svclone/main'

// WGS-only modules (ASCAT + BRASS + HRDetect)
include { ASCAT_ALLELECOUNT           } from '../modules/local/ascat/allelecount/main'
include { ASCAT_RUN                   } from '../modules/local/ascat/run/main'
include { BRASS_GENERATE_BAS          } from '../modules/local/brass/generate_bas/main'
include { BRASS_INPUT                 } from '../modules/local/brass/input/main'
include { BRASS_COVER                 } from '../modules/local/brass/cover/main'
include { BRASS_RUN                   } from '../modules/local/brass/run/main'
include { HRDETECT                    } from '../modules/local/hrdetect/main'
include { SV_SIGNATURES               } from '../modules/local/svsignatures/main'

// Annotation & Signatures modules
include { GERMLINE_ANNOTATE_MAF       } from '../modules/local/germline/annotate_maf/main'
include { SOMATIC_FACETS_ANNOTATION   } from '../modules/local/somatic/facets_annotation/main'
include { GERMLINE_FACETS_ANNOTATION  } from '../modules/local/germline/facets_annotation/main'
include { FACETS_PREVIEW_QC           } from '../modules/local/facets/preview_qc/main'
include { NEOANTIGEN                  } from '../modules/local/neoantigen/main'
include { MUTSIG                      } from '../modules/local/mutsig/main'
include { SPLIT_INTERVALS             } from '../modules/local/splitintervals/main'
include { GERMLINE_COMBINE_HC_VCF     } from '../modules/local/germline/combine_hc_vcf/main'
include { GERMLINE_HARD_FILTER        } from '../modules/local/germline/hard_filter/main'
include { GATK4_MERGEMUTECTSTATS     } from '../modules/local/gatk4/mergemutectstats/main'
include { METADATA_PARSER             } from '../modules/local/metadata_parser/main'

// QC & Reporting modules
include { ALFRED                      } from '../modules/local/alfred/main'
include { CONPAIR_ALL                 } from '../modules/local/conpair/all/main'
include { MULTIQC_SAMPLE              } from '../modules/local/multiqc/sample/main'
include { MULTIQC_SOMATIC             } from '../modules/local/multiqc/somatic/main'
include { MULTIQC_COHORT              } from '../modules/local/multiqc/cohort/main'

// Aggregation modules
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
~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~
*/

workflow TEMPO {

    take:
    ch_input
    ch_fasta
    ch_fasta_fai
    ch_dict
    ch_bwa_index
    ch_dbsnp
    ch_dbsnp_tbi
    ch_known_indels
    ch_known_indels_tbi
    ch_germline_resource
    ch_germline_resource_tbi
    ch_intervals
    ch_pon
    ch_pon_tbi
    ch_bam_input

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ===============
    // WORKFLOW CONTROL FLAGS
    // ===============
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
    // Parse samplesheet and group by patient/sample
    ch_input
        .map { meta, fastq_1, fastq_2 ->
            def new_meta = meta + [ id: meta.sample ]
            [ new_meta, [ fastq_1, fastq_2 ] ]
        }
        .set { ch_reads }

    // ===============
    // RAW READ QC
    // ===============

    FASTQC ( ch_reads )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})

    // ===============
    // PREPROCESSING: TRIM + ALIGN + SORT
    // ===============

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

    // Extract read group ID from FASTQ header
    EXTRACT_READ_GROUP ( FASTP.out.reads )

    FASTP.out.reads
        .join(EXTRACT_READ_GROUP.out.read_group)
        .map { meta, reads, rg_id ->
            def new_meta = meta + [
                read_group: "@RG\\tID:${rg_id.trim()}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA"
            ]
            [ new_meta, reads ]
        }
        .set { ch_reads_with_rg }

    BWAMEM2_MEM (
        ch_reads_with_rg,
        ch_bwa_index,
        ch_fasta,
        true   // sort_bam
    )

    SAMTOOLS_SORT (
        BWAMEM2_MEM.out.bam,
        ch_fasta,
        []
    )

    SAMTOOLS_INDEX_SORTED ( SAMTOOLS_SORT.out.bam )

    // ===============
    // MULTI-LANE MERGE
    // ===============
    // Group BAMs by sample for merging (multi-lane)
    SAMTOOLS_SORT.out.bam
        .map { meta, bam ->
            def new_meta = meta.subMap('patient', 'sample', 'status', 'target') + [id: meta.sample]
            [ new_meta, bam ]
        }
        .groupTuple()
        .branch {
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }
        .set { ch_bams_to_merge }

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

    // ===============
    // MARK DUPLICATES
    // ===============

    GATK4_MARKDUPLICATES (
        ch_bams_merged,
        ch_fasta.map{ it[1] },
        ch_fasta_fai.map{ it[1] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})

    // Index the marked BAMs
    SAMTOOLS_INDEX_MD ( GATK4_MARKDUPLICATES.out.bam )

    // ===============
    // BASE QUALITY SCORE RECALIBRATION
    // ===============

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

    // ===============
    // BAM INPUT (skip alignment for pre-aligned BAMs)
    // ===============
    // BAM inputs bypass FASTQC->FASTP->BWAMEM2->SORT->MERGE->MARKDUP->BQSR
    // and feed directly into the recalibrated BAM channel
    ch_recal_bam_bai
        .mix(ch_bam_input)
        .set { ch_all_recal_bam_bai }

    // ===============
    // TUMOR-NORMAL PAIRING
    // ===============
    ch_all_recal_bam_bai
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
                target:    tumor_meta.target
            ]
            [ pair_meta, tumor_bam, tumor_bai, normal_bam, normal_bai ]
        }
        .set { ch_tumor_normal_pair }

    // ===============
    // SOMATIC SNV/INDEL CALLING
    // ===============

    // ===============
    // MANTA SOMATIC (needed by both SNV and SV workflows)
    // ===============

    if (doWF_manta) {
        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta, nbam, nbai, tbam, tbai, [], [] ]
            }
            .set { ch_manta_input }

        MANTA_SOMATIC (
            ch_manta_input,
            ch_fasta,
            ch_fasta_fai,
            []
        )
    }

    // ===============
    // SPLIT INTERVALS (shared by Mutect2 + HaplotypeCaller scatter-gather)
    // ===============

    if (doWF_SNV || doWF_germSNV) {
        ch_intervals = params.intervals
            ? Channel.value(file(params.intervals, checkIfExists: true))
            : (params.target_intervals
                ? Channel.value(file(params.target_intervals, checkIfExists: true))
                : Channel.value([]))

        SPLIT_INTERVALS (
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_intervals,
            params.scatter_count
        )
        ch_versions = ch_versions.mix(SPLIT_INTERVALS.out.versions)
    }

    // ===============
    // SOMATIC SNV/INDEL CALLING
    // ===============

    if (doWF_SNV) {

        ch_fasta_fai
            .map { meta, fai -> [ meta, fai, [] ] }
            .set { ch_fai_gzi }

        ch_tumor_normal_pair
            .combine(SPLIT_INTERVALS.out.interval_lists.flatten())
            .map { meta, tbam, tbai, nbam, nbai, interval ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}_${interval.baseName}"
                new_meta.original_id = meta.id
                [ new_meta, [tbam, nbam], [tbai, nbai], interval ]
            }
            .set { ch_mutect2_input }

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

        GATK4_MUTECT2.out.vcf
            .map { meta, vcf -> [ meta.original_id ?: meta.id, vcf ] }
            .groupTuple()
            .map { original_id, vcfs ->
                def meta = [id: original_id]
                [ meta, vcfs ]
            }
            .set { ch_mutect2_vcf_gather }

        GATK4_MERGEVCFS (
            ch_mutect2_vcf_gather,
            ch_dict
        )

        GATK4_MUTECT2.out.stats
            .map { meta, stats -> [ meta.original_id ?: meta.id, stats ] }
            .groupTuple()
            .map { original_id, stats ->
                def meta = [id: original_id]
                [ meta, stats ]
            }
            .set { ch_mutect2_stats_gather }

        GATK4_MERGEMUTECTSTATS (
            ch_mutect2_stats_gather
        )

        GATK4_MUTECT2.out.f1r2
            .map { meta, f1r2 -> [ meta.original_id ?: meta.id, f1r2 ] }
            .groupTuple()
            .map { original_id, f1r2s ->
                def meta = [id: original_id]
                [ meta, f1r2s ]
            }
            .set { ch_mutect2_f1r2_gather }

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

        GETPILEUPSUMMARIES_TUMOR.out.table
            .join(GETPILEUPSUMMARIES_NORMAL.out.table)
            .map { meta, tumor_table, normal_table ->
                [ meta, tumor_table, normal_table ]
            }
            .set { ch_contamination_input }

        GATK4_CALCULATECONTAMINATION ( ch_contamination_input )

        GATK4_LEARNREADORIENTATIONMODEL (
            ch_mutect2_f1r2_gather
        )

        GATK4_MERGEVCFS.out.vcf
            .join(GATK4_MERGEVCFS.out.tbi)
            .join(GATK4_MERGEMUTECTSTATS.out.stats)
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
            []
        )

        VCF2MAF (
            ENSEMBLVEP_VEP.out.vcf,
            ch_fasta.map{ it[1] },
            params.vep_cache ? Channel.value(file(params.vep_cache, checkIfExists: true)) : Channel.value([])
        )
        ch_versions = ch_versions.mix(VCF2MAF.out.versions.first())
    }

    // ===============
    // SOMATIC SV CALLING
    // ===============

    if (doWF_SV) {
        // Calls each SV type separately for parallelism and fault tolerance,

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

        DELLY_CALL_SOMATIC.out.vcf
            .groupTuple(by: 0)
            .map { meta, svTypes, vcfs, tbis ->
                [ meta, vcfs, tbis ]
            }
            .set { ch_delly_combine_input }

        DELLY_COMBINE ( ch_delly_combine_input )
        ch_versions = ch_versions.mix(DELLY_COMBINE.out.versions.first())

        SVABA_SOMATIC (
            ch_tumor_normal_pair,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        // Exome: Delly + Manta + SvABA (3 callers)
        // WGS:   + BRASS (4 callers, added in WGS block below)

        if (doWF_manta) {
            DELLY_COMBINE.out.vcf
                .map { meta, vcf, tbi -> [ meta, vcf, tbi, "delly" ] }
                .mix(
                    MANTA_SOMATIC.out.diploid_sv_vcf
                        .join(MANTA_SOMATIC.out.diploid_sv_vcf_tbi)
                        .map { meta, vcf, tbi -> [ meta, vcf, tbi, "manta" ] }
                )
                .mix(
                    SVABA_SOMATIC.out.vcf
                        .map { meta, vcf, tbi -> [ meta, vcf, tbi, "svaba" ] }
                )
                .set { ch_somatic_sv_callers_base }
        }
    }

    // ===============
    // FACETS (Copy Number)
    // ===============

    if (doWF_facets) {
        SNPPILEUP (
            ch_tumor_normal_pair,
            params.facets_vcf ? Channel.value(file(params.facets_vcf, checkIfExists: true)) : Channel.value([])
        )
        ch_versions = ch_versions.mix(SNPPILEUP.out.versions.first())

        FACETS ( SNPPILEUP.out.pileup )
        ch_versions = ch_versions.mix(FACETS.out.versions.first())
    }

    // ===============
    // FACETS PREVIEW QC & ANNOTATION
    // ===============

    if (doWF_facets) {
        FACETS_PREVIEW_QC (
            FACETS.out.facets_output
                .join(SNPPILEUP.out.pileup)
                .map { meta, output, pileup -> [ meta, output, pileup ] }
        )

        if (doWF_SNV) {
            FACETS.out.hisens_rdata
                .join(VCF2MAF.out.maf)
                .set { ch_somatic_facets_anno_input }

            SOMATIC_FACETS_ANNOTATION (
                ch_somatic_facets_anno_input
            )
        }

    }

    // ===============
    // MUTATION SIGNATURES
    // ===============

    if (doWF_mutSig) {
        if (doWF_facets && doWF_SNV) {
            MUTSIG ( SOMATIC_FACETS_ANNOTATION.out.final_maf )
        }
    }

    // ===============
    // MSI
    // ===============

    if (doWF_msiSensor) {
        MSISENSORPRO_SCAN ( ch_fasta )

        //              path(tumor), path(tumor_index), path(intervals)

        //        path(msisensor_scan)

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

    // ===============
    // HLA TYPING & LOH
    // ===============

    if (doWF_loh) {
        POLYSOLVER (
            ch_recal_branched.normal.map { meta, bam, bai -> [ meta, bam, bai ] }
        )
        ch_versions = ch_versions.mix(POLYSOLVER.out.versions.first())

        //              path(normal_bam), path(normal_bai), path(hla_types)

        ch_tumor_normal_pair
            .map { meta, tbam, tbai, nbam, nbai ->
                [ meta.normal_id, meta, tbam, tbai, nbam, nbai ]
            }
            .combine(
                POLYSOLVER.out.hla_types.map { meta, hla -> [ meta.sample ?: meta.id, hla ] },
                by: 0
            )
            .map { normal_id, meta, tbam, tbai, nbam, nbai, hla_types ->
                [ meta.id, meta, tbam, tbai, nbam, nbai, hla_types ]
            }
            .join(
                FACETS.out.purity.map { meta, out_file -> [ meta.id, out_file ] }
            )
            .map { id, meta, tbam, tbai, nbam, nbai, hla_types, purity_out ->
                [ meta, tbam, tbai, nbam, nbai, hla_types, purity_out ]
            }
            .set { ch_lohhla_input }

        LOHHLA (
            ch_lohhla_input,
            params.hla_fasta ? Channel.value(file(params.hla_fasta, checkIfExists: true)) : Channel.value([]),
            params.hla_dat   ? Channel.value(file(params.hla_dat, checkIfExists: true))   : Channel.value([])
        )
        ch_versions = ch_versions.mix(LOHHLA.out.versions.first())
    }

    // ===============
    // GERMLINE VARIANT CALLING
    // ===============

    if (doWF_germSNV) {
        ch_recal_branched.normal
            .combine(SPLIT_INTERVALS.out.interval_lists.flatten())
            .map { meta, bam, bai, interval ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}_${interval.baseName}"
                new_meta.original_id = meta.id
                [ new_meta, bam, bai, interval, [] ]
            }
            .set { ch_hc_input }

        GATK4_HAPLOTYPECALLER (
            ch_hc_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_tbi
        )

        // Groups VCFs by original sample ID, then concat + normalize + dedup

        GATK4_HAPLOTYPECALLER.out.vcf
            .map { meta, vcf -> [ meta.original_id ?: meta.id, vcf ] }
            .groupTuple()
            .join(
                GATK4_HAPLOTYPECALLER.out.tbi
                    .map { meta, tbi -> [ meta.original_id ?: meta.id, tbi ] }
                    .groupTuple()
            )
            .map { original_id, vcfs, tbis ->
                def meta = [id: original_id, sample: original_id]
                [ meta, vcfs, tbis ]
            }
            .set { ch_hc_combine_input }

        GERMLINE_COMBINE_HC_VCF (
            ch_hc_combine_input,
            ch_fasta,
            ch_fasta_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(GERMLINE_COMBINE_HC_VCF.out.versions.first())

        GERMLINE_HARD_FILTER (
            GERMLINE_COMBINE_HC_VCF.out.vcf,
            ch_fasta,
            ch_fasta_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(GERMLINE_HARD_FILTER.out.versions.first())

        //        path(call_regions), path(call_regions_tbi)

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

        GERMLINE_HARD_FILTER.out.vcf
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

    // ===============
    // GERMLINE SV CALLING
    // ===============

    if (doWF_germSV) {
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

        DELLY_CALL_GERMLINE.out.vcf
            .groupTuple(by: 0)
            .map { meta, svTypes, vcfs, tbis -> [ meta, vcfs, tbis ] }
            .set { ch_delly_combine_germline_input }

        DELLY_COMBINE_GERMLINE ( ch_delly_combine_germline_input )
        ch_versions = ch_versions.mix(DELLY_COMBINE_GERMLINE.out.versions.first())

        SVABA_GERMLINE (
            ch_recal_branched.normal,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        DELLY_COMBINE_GERMLINE.out.vcf
            .map { meta, vcf, tbi -> [ meta, vcf, tbi, "delly" ] }
            .mix(
                MANTA_GERMLINE.out.vcf
                    .map { meta, vcf, tbi -> [ meta, vcf, tbi, "manta" ] }
            )
            .mix(
                SVABA_GERMLINE.out.vcf
                    .map { meta, vcf, tbi -> [ meta, vcf, tbi, "svaba" ] }
            )
            .groupTuple(by: 0, size: 3)
            .map { meta, vcfs, tbis, callers ->
                [ meta, vcfs, tbis, callers ]
            }
            .set { ch_germline_sv_merge_input }

        GERMLINE_MERGE_SV ( ch_germline_sv_merge_input )
        ch_versions = ch_versions.mix(GERMLINE_MERGE_SV.out.versions.first())

        SVTOOLS_VCF2BEDPE_GERMLINE (
            GERMLINE_MERGE_SV.out.vcf
        )

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

    // ===============
    // GERMLINE FACETS ANNOTATION (after germline SNV so GERMLINE_COMBINE_CHANNEL is available)
    // ===============

    if (doWF_facets && doWF_germSNV) {
        FACETS.out.hisens_rdata
            .join(GERMLINE_ANNOTATE_MAF.out.maf_file)
            .set { ch_germline_facets_anno_input }

        GERMLINE_FACETS_ANNOTATION (
            ch_germline_facets_anno_input
        )
    }

    // ===============
    // NEOANTIGEN PREDICTION (after LOH/POLYSOLVER so POLYSOLVER.out is available)
    // ===============

    if (doWF_loh && doWF_SNV) {
        ch_neoantigen_cdna = params.neoantigen_cdna
            ? Channel.value(file(params.neoantigen_cdna, checkIfExists: true))
            : Channel.value([])
        ch_neoantigen_cds = params.neoantigen_cds
            ? Channel.value(file(params.neoantigen_cds, checkIfExists: true))
            : Channel.value([])

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

    // ===============
    // METADATA PARSER (after MSI + LOH + MUTSIG so all inputs available)
    // ===============

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

    // ===============
    // WGS-ONLY: ASCAT + BRASS + HRDetect
    // ===============

    if (isWGS && doWF_SV) {
        ASCAT_ALLELECOUNT (
            ch_tumor_normal_pair,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.snp_gc_corrections ? Channel.value(file(params.snp_gc_corrections, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

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

        // Run separately on tumor and normal BAMs

        ch_all_recal_bam_bai
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

        BRASS_COVER (
            ch_brass_input_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            params.brass_ref_dir ? Channel.value(file(params.brass_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE')),
            params.vagrent_ref_dir ? Channel.value(file(params.vagrent_ref_dir, checkIfExists: true)) : Channel.value(file('NO_FILE'))
        )

        // Combine BRASS with BAMs + ASCAT
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

        // WGS: Add BRASS to the SV callers channel for 4-caller merge

        if (doWF_manta) {
            ch_somatic_sv_callers_base
                .mix(
                    BRASS_RUN.out.brass_vcf
                        .map { meta, vcf, tbi -> [ meta, vcf, tbi, "brass" ] }
                )
                .set { ch_somatic_sv_callers_with_brass }
        }
    }

    // ===============
    // SOMATIC SV MERGE + ANNOTATION
    // Runs after WGS block so BRASS is available for WGS
    // Exome: 3 callers (Delly+Manta+SvABA), WGS: 4 callers (+BRASS)
    // ===============

    if (doWF_SV && doWF_manta) {
        def sv_caller_count = isWGS ? 4 : 3
        def ch_sv_callers_final = isWGS ? ch_somatic_sv_callers_with_brass : ch_somatic_sv_callers_base

        ch_sv_callers_final
            .groupTuple(by: 0, size: sv_caller_count)
            .map { meta, vcfs, tbis, callers -> [ meta, vcfs, tbis, callers ] }
            .set { ch_sv_merge_input }

        SOMATIC_MERGE_SV ( ch_sv_merge_input )
        ch_versions = ch_versions.mix(SOMATIC_MERGE_SV.out.versions.first())

        // VCF2BEDPE somatic
        SVTOOLS_VCF2BEDPE_SOMATIC ( SOMATIC_MERGE_SV.out.vcf )

        // iAnnotateSV somatic
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
            ch_repeat_masker, ch_mapability_blacklist,
            ch_sv_blacklist_bed, ch_sv_blacklist_bedpe,
            ch_sv_blacklist_foldback_bedpe, ch_sv_blacklist_te_bedpe,
            ch_splice_sites, params.genome ?: 'GRCh37'
        )

        // ClusterSV (cluster breakpoints)
        CLUSTERSV ( IANNOTATESV_SOMATIC.out.bedpe_pass, params.genome ?: 'GRCh37' )

        // WGS-only: HRDetect + SV Signatures
        if (isWGS) {
            if (doWF_SNV && doWF_facets) {
                ch_hrdetect_script = params.hrdetect_script
                    ? Channel.value(file(params.hrdetect_script, checkIfExists: true))
                    : Channel.value(file('NO_FILE'))
                SOMATIC_FACETS_ANNOTATION.out.final_maf
                    .join(FACETS.out.hisens_seg)
                    .join(IANNOTATESV_SOMATIC.out.bedpe_pass)
                    .map { meta, maf, cnv, bedpe -> [ meta, maf, cnv, bedpe ] }
                    .set { ch_hrdetect_input }
                HRDETECT ( ch_hrdetect_input, ch_hrdetect_script )
            }

            ch_svsig_script = params.sv_signature_script
                ? Channel.value(file(params.sv_signature_script, checkIfExists: true))
                : Channel.value(file('NO_FILE'))
            SV_SIGNATURES ( IANNOTATESV_SOMATIC.out.bedpe_pass, ch_svsig_script )
        }
    }

    // ===============
    // SVCircos (circos plot visualization)
    // ===============

    if (doWF_SV && doWF_facets && doWF_manta) {
        IANNOTATESV_SOMATIC.out.bedpe_pass
            .join(FACETS.out.hisens_seg)
            .set { ch_svcircos_input }
        SVCIRCOS ( ch_svcircos_input, params.genome ?: 'GRCh37' )
    }

    // ===============
    // SVCLONE
    // ===============

    if (doWF_SV && doWF_facets && doWF_SNV && doWF_manta) {
        ch_tumor_normal_pair
            .join(IANNOTATESV_SOMATIC.out.bedpe_pass)
            .join(SOMATIC_FACETS_ANNOTATION.out.final_maf)
            .join(FACETS.out.hisens_seg)
            .join(FACETS.out.purity)
            .set { ch_svclone_input }
        SVCLONE ( ch_svclone_input )
    }

    // ===============
    // QC
    // ===============

    if (doWF_QC) {
        ch_all_recal_bam_bai
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

        QUALIMAP_BAMQC (
            ch_all_recal_bam_bai.map { meta, bam, bai -> [ meta, bam ] },
            []
        )

        CONPAIR_PILEUP (
            ch_all_recal_bam_bai.map { meta, bam, bai -> [ meta, bam, bai ] },
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )
        ch_versions = ch_versions.mix(CONPAIR_PILEUP.out.versions.first())

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

        CONPAIR_ALL (
            ch_conpair_concordance_input,
            ch_fasta.map{ it[1] },
            ch_fasta_fai.map{ it[1] },
            ch_dict.map{ it[1] }
        )

        ch_all_recal_bam_bai
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

        ch_multiqc_sample_configs = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

        def ch_fastp_json_per_sample = FASTP.out.json
            .map { meta, json -> [ meta.sample, json ] }
            .groupTuple()
        ALFRED.out.alfred_qc
            .map { meta, rg_n, rg_y ->
                [ meta.sample, meta, rg_n, rg_y ]
            }
            .join( ch_fastp_json_per_sample, by: 0 )
            .map { sample, meta, rg_n, rg_y, fastp_jsons ->
                [ meta, rg_n, rg_y, fastp_jsons ]
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

        // Per tumor-normal pair QC metrics aggregation

        if (doWF_SNV && doWF_facets) {
            def ch_qualimap_by_sample = QUALIMAP_BAMQC.out.results
                .map { meta, dir -> [ meta.sample, dir ] }
            CONPAIR_CONCORDANCE.out.concordance
                .map { meta, conc -> [ meta.tumor_id, meta, conc ] }
                .combine( ch_qualimap_by_sample, by: 0 )
                .map { tumor_id, meta, conc, qualimap_tumor ->
                    [ meta.normal_id, meta, conc, qualimap_tumor ]
                }
                .combine( ch_qualimap_by_sample, by: 0 )
                .map { normal_id, meta, conc, qualimap_tumor, qualimap_normal ->
                    [ meta, conc, qualimap_tumor, qualimap_normal ]
                }
                .join(FACETS.out.summary_out)
                .join(FACETS_PREVIEW_QC.out.facets_qc)
                .map { meta, conpair, qualimap_t, qualimap_n, facets_sum, facets_qc ->
                    [ meta, conpair, qualimap_t, qualimap_n, facets_sum, facets_qc ]
                }
                .set { ch_multiqc_somatic_input }

            MULTIQC_SOMATIC (
                ch_multiqc_somatic_input,
                ch_multiqc_sample_configs.toList()
            )
        }
    }

    // ===============
    // MultiQC (nf-core - fallback aggregation)
    // ===============

    if (doWF_QC && !params.skip_multiqc) {

        ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
        ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)   : Channel.empty()

        //        path(multiqc_config)
        //        path(extra_multiqc_config)
        //        path(multiqc_logo)
        //        path(replace_names)
        //        path(sample_names)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],  // replace_names
            []
        )
    }

    // ===============
    // AGGREGATION (cohort-level outputs)
    // ===============
    if (params.aggregate) {
        def aggregateIsFile = params.aggregate instanceof String && params.aggregate != 'true' && file(params.aggregate).exists()
        if (aggregateIsFile) {
            Channel.fromPath(params.aggregate)
                .splitCsv(sep: '\t', header: true)
                .map { row -> [ row.COHORT, row.TUMOR_ID, row.NORMAL_ID, row.PATH ?: '' ] }
                .set { ch_aggregate_raw }

            ch_aggregate_raw
                .map { cohort, tid, nid, path -> [ cohort, tid, nid ] }
                .set { ch_aggregate_map }

            // Rows with PATH column — used for aggregate-only file resolution
            ch_aggregate_raw
                .filter { cohort, tid, nid, path -> path }
                .set { ch_aggregate_with_path }
        } else {
            ch_tumor_normal_pair
                .map { meta, tbam, tbai, nbam, nbai -> [ "default_cohort", meta.tumor_id, meta.normal_id ] }
                .set { ch_aggregate_map }
            ch_aggregate_with_path = Channel.empty()
        }
        if (doWF_SNV && doWF_facets) {
            def ch_maf_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/somatic/${tid}__${nid}/*/*.final.maf")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ tid, nid, f ]
                }
            def ch_maf_by_pair = SOMATIC_FACETS_ANNOTATION.out.final_maf
                .map { meta, maf -> [ meta.tumor_id, meta.normal_id, maf ] }
                .mix(ch_maf_from_path)
            ch_aggregate_map.combine(ch_maf_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, mafs -> [ cohort, mafs ] }
                .set { ch_agg_somatic_maf }
            ch_agg_somatic_maf.map { cohort, mafs -> cohort }.set { ch_agg_maf_cohort }
            ch_agg_somatic_maf.map { cohort, mafs -> mafs }.set { ch_agg_maf_files }
            AGGREGATE_SOMATIC_MAF ( ch_agg_maf_cohort, ch_agg_maf_files )
        }
        if (doWF_SV && doWF_manta) {
            def ch_sv_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/somatic/${tid}__${nid}/*/*.delly.manta.vcf.gz")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ tid, nid, f ]
                }
            def ch_sv_by_pair = IANNOTATESV_SOMATIC.out.bedpe_pass
                .map { meta, bedpe -> [ meta.tumor_id, meta.normal_id, bedpe ] }
                .mix(ch_sv_from_path)
            ch_aggregate_map.combine(ch_sv_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_somatic_sv }
            ch_agg_somatic_sv.map { it[0] }.set { ch_agg_sv_cohort }
            ch_agg_somatic_sv.map { it[1] }.set { ch_agg_sv_files }
            AGGREGATE_SOMATIC_SV ( ch_agg_sv_cohort, ch_agg_sv_files )
        }
        if (doWF_facets) {
            def ch_facets_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def pur = file("${path}/somatic/${tid}__${nid}/*/*/*/*_purity.seg")
                    def his = file("${path}/somatic/${tid}__${nid}/*/*/*/*_hisens.seg")
                    def out_f = file("${path}/somatic/${tid}__${nid}/*/*/*_OUT.txt")
                    def arm = file("${path}/somatic/${tid}__${nid}/*/*/*/*.arm_level.txt")
                    def gene = file("${path}/somatic/${tid}__${nid}/*/*/*/*.gene_level.txt")
                    [ tid, nid,
                      pur instanceof List ? (pur.size() > 0 ? pur[0] : file('NO_FILE')) : pur,
                      his instanceof List ? (his.size() > 0 ? his[0] : file('NO_FILE')) : his,
                      out_f instanceof List ? (out_f.size() > 0 ? out_f[0] : file('NO_FILE')) : out_f,
                      arm instanceof List ? (arm.size() > 0 ? arm[0] : file('NO_FILE')) : arm,
                      gene instanceof List ? (gene.size() > 0 ? gene[0] : file('NO_FILE')) : gene ]
                }
            def ch_facets_by_pair = FACETS.out.purity_seg
                .join(FACETS.out.hisens_seg).join(FACETS.out.summary_out)
                .join(FACETS.out.arm_level).join(FACETS.out.gene_level)
                .map { meta, pur, his, out, arm, gene -> [ meta.tumor_id, meta.normal_id, pur, his, out, arm, gene ] }
                .mix(ch_facets_from_path)
            ch_aggregate_map.combine(ch_facets_by_pair, by: [1,2])
                .groupTuple(by: 0)
                .map { cohort, tids, nids, purs, hiss, outs, arms, genes -> [ cohort, purs, hiss, outs, arms, genes ] }
                .set { ch_agg_facets }
            AGGREGATE_SOMATIC_FACETS ( ch_agg_facets )
        }
        if (doWF_loh && doWF_SNV) {
            def ch_neo_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/somatic/${tid}__${nid}/*/*.all_neoantigen_predictions.txt")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ tid, nid, f ]
                }
            def ch_neo_by_pair = NEOANTIGEN.out.predictions
                .map { meta, neo -> [ meta.tumor_id, meta.normal_id, neo ] }
                .mix(ch_neo_from_path)
            ch_aggregate_map.combine(ch_neo_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_netmhc }
            ch_agg_netmhc.map { it[0] }.set { ch_agg_netmhc_cohort }
            ch_agg_netmhc.map { it[1] }.set { ch_agg_netmhc_files }
            AGGREGATE_SOMATIC_NETMHC ( ch_agg_netmhc_cohort, ch_agg_netmhc_files )
        }
        if (doWF_mdParse) {
            def ch_md_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/somatic/${tid}__${nid}/*/*.sample_data.txt")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ tid, nid, f ]
                }
            def ch_md_by_pair = METADATA_PARSER.out.metadata
                .map { meta, md -> [ meta.tumor_id, meta.normal_id, md ] }
                .mix(ch_md_from_path)
            ch_aggregate_map.combine(ch_md_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_metadata }
            ch_agg_metadata.map { it[0] }.set { ch_agg_md_cohort }
            ch_agg_metadata.map { it[1] }.set { ch_agg_md_files }
            AGGREGATE_SOMATIC_METADATA ( ch_agg_md_cohort, ch_agg_md_files )
        }
        if (doWF_loh) {
            def ch_lohhla_pred_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def pred = file("${path}/somatic/${tid}__${nid}/*/*.DNA.HLAlossPrediction_CI.txt")
                    def cpn = file("${path}/somatic/${tid}__${nid}/*/*DNA.IntegerCPN_CI.txt")
                    def p = pred instanceof List ? (pred.size() > 0 ? pred[0] : file('NO_FILE')) : pred
                    def c = cpn instanceof List ? (cpn.size() > 0 ? cpn[0] : file('NO_FILE')) : cpn
                    [ tid, nid, p, c ]
                }
            def ch_lohhla_by_pair = LOHHLA.out.predictions
                .join(LOHHLA.out.integer_cpn)
                .map { meta, pred, cpn -> [ meta.tumor_id, meta.normal_id, pred, cpn ] }
                .mix(ch_lohhla_pred_from_path)
            ch_aggregate_map.combine(ch_lohhla_by_pair, by: [1,2])
                .groupTuple(by: 0)
                .map { cohort, tids, nids, preds, cpns -> [ cohort, preds, cpns ] }
                .set { ch_agg_lohhla }
            ch_agg_lohhla.map { it[0] }.set { ch_agg_lohhla_cohort }
            ch_agg_lohhla.map { it[1] }.set { ch_agg_lohhla_preds }
            ch_agg_lohhla.map { it[2] }.set { ch_agg_lohhla_cpns }
            AGGREGATE_SOMATIC_LOHHLA ( ch_agg_lohhla_cohort, ch_agg_lohhla_preds, ch_agg_lohhla_cpns )
        }
        if (isWGS && doWF_SV && doWF_SNV && doWF_facets) {
            def ch_hrd_by_pair = HRDETECT.out.hrdetect_output
                .map { meta, hrd -> [ meta.tumor_id, meta.normal_id, hrd ] }
            ch_aggregate_map.combine(ch_hrd_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_hrdetect }
            ch_agg_hrdetect.map { it[0] }.set { ch_agg_hrd_cohort }
            ch_agg_hrdetect.map { it[1] }.set { ch_agg_hrd_files }
            AGGREGATE_SOMATIC_HRDETECT ( ch_agg_hrd_cohort, ch_agg_hrd_files )
        }
        if (doWF_SV && doWF_facets && doWF_SNV && doWF_manta) {
            def ch_svc_by_pair = SVCLONE.out.cluster_certainty
                .map { meta, sv_cert, snv_cert -> [ meta.tumor_id, meta.normal_id, sv_cert, snv_cert ] }
            ch_aggregate_map.combine(ch_svc_by_pair, by: [1,2])
                .groupTuple(by: 0)
                .map { cohort, tids, nids, svs, snvs -> [ cohort, svs, snvs ] }
                .set { ch_agg_svclone }
            ch_agg_svclone.map { it[0] }.set { ch_agg_svc_cohort }
            ch_agg_svclone.map { it[1] }.set { ch_agg_svc_sv }
            ch_agg_svclone.map { it[2] }.set { ch_agg_svc_snv }
            AGGREGATE_SOMATIC_SVCLONE ( ch_agg_svc_cohort, ch_agg_svc_sv, ch_agg_svc_snv )
        }
        if (isWGS && doWF_SV && doWF_manta) {
            def ch_svsig_by_pair = SV_SIGNATURES.out.catalogues
                .join(SV_SIGNATURES.out.sv_signatures)
                .map { meta, pdf, sigs -> [ meta.tumor_id, meta.normal_id, pdf, sigs ] }
            ch_aggregate_map.combine(ch_svsig_by_pair, by: [1,2])
                .groupTuple(by: 0)
                .map { cohort, tids, nids, pdfs, sigs -> [ cohort, pdfs, sigs ] }
                .set { ch_agg_svsig }
            ch_agg_svsig.map { it[0] }.set { ch_agg_svsig_cohort }
            ch_agg_svsig.map { it[1] }.set { ch_agg_svsig_pdfs }
            ch_agg_svsig.map { it[2] }.set { ch_agg_svsig_exps }
            AGGREGATE_SOMATIC_SVSIGNATURES ( ch_agg_svsig_cohort, ch_agg_svsig_pdfs, ch_agg_svsig_exps )
        }
        if (doWF_germSNV) {
            def ch_gmaf_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/germline/${nid}/*/${tid}__${nid}.germline.final.maf")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ tid, nid, f ]
                }
            def ch_gmaf_by_pair = GERMLINE_ANNOTATE_MAF.out.maf_file
                .map { meta, maf -> [ meta.tumor_id, meta.normal_id, maf ] }
                .mix(ch_gmaf_from_path)
            ch_aggregate_map.combine(ch_gmaf_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_germline_maf }
            ch_agg_germline_maf.map { it[0] }.set { ch_agg_gmaf_cohort }
            ch_agg_germline_maf.map { it[1] }.set { ch_agg_gmaf_files }
            AGGREGATE_GERMLINE_MAF ( ch_agg_gmaf_cohort, ch_agg_gmaf_files )
        }
        if (doWF_germSV) {
            def ch_gsv_from_path = ch_aggregate_with_path
                .map { cohort, tid, nid, path ->
                    def resolved = file("${path}/germline/${nid}/*/*.delly.manta.vcf.gz")
                    def f = resolved instanceof List ? (resolved.size() > 0 ? resolved[0] : file('NO_FILE')) : resolved
                    [ nid, cohort, f ]
                }
            def ch_gsv_by_normal = IANNOTATESV_GERMLINE.out.bedpe_pass
                .map { meta, bedpe -> [ meta.sample ?: meta.id, bedpe ] }
            def ch_gsv_from_pipeline = ch_aggregate_map.map { cohort, tid, nid -> [ nid, cohort ] }
                .combine(ch_gsv_by_normal, by: 0)
                .map { nid, cohort, bedpe -> [ cohort, bedpe ] }
            ch_gsv_from_pipeline
                .mix(ch_gsv_from_path.map { nid, cohort, f -> [ cohort, f ] })
                .groupTuple(by: 0).map { cohort, files -> [ cohort, files.unique() ] }
                .set { ch_agg_germline_sv }
            ch_agg_germline_sv.map { it[0] }.set { ch_agg_gsv_cohort }
            ch_agg_germline_sv.map { it[1] }.set { ch_agg_gsv_files }
            AGGREGATE_GERMLINE_SV ( ch_agg_gsv_cohort, ch_agg_gsv_files )
        }
        if (doWF_QC) {
            ALFRED.out.alfred_qc
                .map { meta, rg_n, rg_y -> [ rg_n, rg_y ] }
                .flatten().collect()
                .set { ch_agg_alfred }
            PICARD_COLLECTHSMETRICS.out.metrics
                .map { meta, m -> m }.collect()
                .set { ch_agg_hsmetrics }
            AGGREGATE_QC_BAM ( Channel.value("default_cohort"), ch_agg_alfred, ch_agg_hsmetrics )
        }
        if (doWF_QC) {
            def ch_conc_by_pair = CONPAIR_CONCORDANCE.out.concordance
                .map { meta, c -> [ meta.tumor_id, meta.normal_id, c ] }
            def ch_cont_by_pair = CONPAIR_CONCORDANCE.out.contamination
                .map { meta, c -> [ meta.tumor_id, meta.normal_id, c ] }
            ch_aggregate_map.combine(ch_conc_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_conc }
            ch_aggregate_map.combine(ch_cont_by_pair, by: [1,2])
                .groupTuple(by: 0).map { cohort, tids, nids, files -> [ cohort, files ] }
                .set { ch_agg_cont }
            ch_agg_conc.join(ch_agg_cont, by: 0)
                .set { ch_agg_conpair }
            ch_agg_conpair.map { it[0] }.set { ch_agg_conpair_cohort }
            ch_agg_conpair.map { it[1] }.set { ch_agg_conpair_conc }
            ch_agg_conpair.map { it[2] }.set { ch_agg_conpair_cont }
            AGGREGATE_QC_CONPAIR ( ch_agg_conpair_cohort, ch_agg_conpair_conc, ch_agg_conpair_cont )
        }
        MULTIQC_SAMPLE.out.multiqc_report
            .map { meta, html, data -> [ html, data ] }.flatMap()
            .mix(doWF_QC ? CONPAIR_ALL.out.conpair_output.map { meta, conc, cont -> [ conc, cont ] }.flatMap() : Channel.empty())
            .collect().set { ch_cohort_multiqc_input }
        MULTIQC_COHORT ( ch_cohort_multiqc_input )
        ch_versions = ch_versions.mix(MULTIQC_COHORT.out.versions)
    }

    emit:
    versions       = ch_versions
    multiqc_report = MULTIQC.out.report
}

/*
~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~
*/
