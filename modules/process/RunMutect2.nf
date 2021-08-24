process RunMutect2 {
  tag {idTumor + "__" + idNormal + "@" + intervalBed.baseName}

  input:
    tuple val(id), val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(intervalBed)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)
    val(tools)
    val(runSomatic)

  output:
    tuple val(id), val(idTumor), val(idNormal), val(target), path("*filtered.vcf.gz"), path("*filtered.vcf.gz.tbi"), path("*Mutect2FilteringStats.tsv"), emit: forMutect2Combine

  when: "mutect2" in tools && runSomatic

  script:
  mutect2Vcf = "${idTumor}__${idNormal}_${intervalBed.baseName}.vcf.gz"
  prefix = "${mutect2Vcf}".replaceFirst(".vcf.gz", "")
  """
  gatk --java-options -Xmx8g \
    Mutect2 \
    --reference ${genomeFile} \
    --intervals ${intervalBed} \
    --input ${bamTumor} \
    --tumor-sample ${idTumor} \
    --input ${bamNormal} \
    --normal-sample ${idNormal} \
    --output ${mutect2Vcf}

  gatk --java-options -Xmx8g \
    FilterMutectCalls \
    --variant ${mutect2Vcf} \
    --stats ${prefix}.Mutect2FilteringStats.tsv \
    --output ${prefix}.filtered.vcf.gz
  """
}
