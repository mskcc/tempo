process SomaticAnnotateMaf {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/combined_mutations/", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfMerged)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict), path(vepCache), path(isoforms)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.maf"), emit: mafFile
    path("${outputPrefix}.unfiltered.maf"), emit: unfilteredMafOutput

  script:
  outputPrefix = "${idTumor}__${idNormal}.somatic"
  mutect2InfoCols = "MBQ,MFRL,MMQ,MPOS,OCM,RPA,STR,ECNT"
  strelka2InfoCols = "RU,IC,MQ,SNVSB"
  strelka2FormatCols = "FDP,SUBDP"
  formatCols = "alt_count_raw,alt_count_raw_fwd,alt_count_raw_rev,ref_count_raw,ref_count_raw_fwd,ref_count_raw_rev,depth_raw,depth_raw_fwd,depth_raw_rev"
  formatCols = formatCols + "," + strelka2FormatCols
  if (target == "wgs") {
    infoCols = "MuTect2,Strelka2,Custom_filters,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,AC,AF,AC_nfe_seu,AF_nfe_seu,AC_afr,AF_afr,AC_nfe_onf,AF_nfe_onf,AC_amr,AF_amr,AC_eas,AF_eas,AC_nfe_nwe,AF_nfe_nwe,AC_nfe_est,AF_nfe_est,AC_nfe,AF_nfe,AC_fin,AF_fin,AC_asj,AF_asj,AC_oth,AF_oth,AC_popmax,AN_popmax,AF_popmax"
    infoCols = infoCols + "," + mutect2InfoCols + "," + strelka2InfoCols
  }
  else {
    infoCols = "MuTect2,Strelka2,Custom_filters,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax"
    infoCols = infoCols + "," + mutect2InfoCols + "," + strelka2InfoCols
  }
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 10 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --retain-fmt ${formatCols} \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0
  
  python /usr/bin/oncokb_annotator/MafAnnotator.py \
    -u "https://data-legacy.oncokb.aws.mskcc.org/api/v1" \
    -i ${outputPrefix}.raw.maf \
    -o ${outputPrefix}.raw.oncokb.maf

  Rscript --no-init-file /usr/bin/filter-somatic-maf.R \
    --tumor-vaf ${params.somaticVariant.tumorVaf} \
    --tumor-depth ${params.somaticVariant.tumorDepth} \
    --tumor-count ${params.somaticVariant.tumorCount} \
    --normal-depth ${params.somaticVariant.normalDepth} \
    --normal-count ${params.somaticVariant.normalCount} \
    --gnomad-allele-frequency ${params.somaticVariant.gnomadAf} \
    --normal-panel-count ${params.somaticVariant.ponCount} \
    --maf-file ${outputPrefix}.raw.oncokb.maf \
    --output-prefix ${outputPrefix}
  """
}
