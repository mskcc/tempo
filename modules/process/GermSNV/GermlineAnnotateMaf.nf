process GermlineAnnotateMaf {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/combined_mutations", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfMerged)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict), path(vepCache), path(isoforms)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.maf"), emit: mafFileGermline
    path("${outputPrefix}.unfiltered.maf"), emit: unfilteredMafOutputGermline

  script:
  outputPrefix = "${idTumor}__${idNormal}.germline"
  if (target == 'wgs') {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,AC,AF,AC_nfe_seu,AF_nfe_seu,AC_afr,AF_afr,AC_nfe_onf,AF_nfe_onf,AC_amr,AF_amr,AC_eas,AF_eas,AC_nfe_nwe,AF_nfe_nwe,AC_nfe_est,AF_nfe_est,AC_nfe,AF_nfe,AC_fin,AF_fin,AC_asj,AF_asj,AC_oth,AF_oth,AC_popmax,AN_popmax,AF_popmax"
  }
  else {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax"
  }
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 4 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0

  Rscript --no-init-file /usr/bin/filter-germline-maf.R \
    --normal-depth ${params.germlineVariant.normalDepth} \
    --normal-vaf ${params.germlineVariant.normalVaf} \
    --maf-file ${outputPrefix}.raw.maf \
    --output-prefix ${outputPrefix}
  """
}
