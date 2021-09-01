params.outDir = ""

process SomaticDellyCall {
  tag {idTumor + "__" + idNormal + '@' + svType}

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  
  input:
    each svType
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex), path(svCallingExcludeRegions)
    val(tools)
    val(runSomatic)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi"), emit: dellyFilter4Combine
    tuple path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi"), emit: dellyOutput

  when: "delly" in tools && runSomatic

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idTumor}__${idNormal}_${svType}.bcf \
    ${bamTumor} ${bamNormal}

  echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv

  delly filter \
    --filter somatic \
    --samples samples.tsv \
    --outfile ${idTumor}__${idNormal}_${svType}.filter.bcf \
    ${idTumor}__${idNormal}_${svType}.bcf


  bcftools view --output-type z ${idTumor}__${idNormal}_${svType}.filter.bcf > ${idTumor}__${idNormal}_${svType}.delly.vcf.gz
  tabix --preset vcf ${idTumor}__${idNormal}_${svType}.delly.vcf.gz
  """
}
