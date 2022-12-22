process SomaticDellyCall {
  tag "${idTumor + "__" + idNormal + '@' + svType}"

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  
  input:
    each svType
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex), path(svCallingExcludeRegions)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi"), emit: dellyFilter4Combine
    tuple path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), path("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi"), emit: dellyOutput

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
    -a .05 \
    --samples samples.tsv \
    --outfile ${idTumor}__${idNormal}_${svType}.filter.bcf \
    ${idTumor}__${idNormal}_${svType}.bcf

  # Filter variants that have low supporting reads in the tumor or high supporting reads in the normal
  # DV = discordant reads
  # RV = split reads
  bcftools view \\
    -s ${idTumor},${idNormal} \\
    ${idTumor}__${idNormal}_${svType}.filter.bcf | \\
  bcftools filter \\
    --soft-filter tumor_read_supp -m + \\
    -e "FORMAT/DV[0] < 5 | FORMAT/RV[0] < 2" | \\
  bcftools filter \\
    --soft-filter normal_read_supp -m + \\
    -e "FORMAT/DV[1] > 0 | FORMAT/RV[1] > 0" | \\
  bcftools view --output-type z > \\
    ${idTumor}__${idNormal}_${svType}.delly.vcf.gz

  tabix --preset vcf ${idTumor}__${idNormal}_${svType}.delly.vcf.gz
  """
}
