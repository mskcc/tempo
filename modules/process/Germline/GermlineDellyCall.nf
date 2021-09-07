process GermlineDellyCall {
  tag {idNormal + '@' + svType}

  publishDir "${params.outDir}/germline/${idNormal}/delly", mode: params.publishDirMode, pattern: "*delly.vcf.{gz,gz.tbi}"

  input:
    each svType
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex), path(svCallingExcludeRegions)
    val(tools)
    val(runGermline)

  output:
    tuple val(idNormal), val(target), path("${idNormal}_${svType}.delly.vcf.gz"), path("${idNormal}_${svType}.delly.vcf.gz.tbi"), emit: dellyFilter4CombineGermline
    tuple path("*delly.vcf.gz"), path("*delly.vcf.gz.tbi"), emit: dellyOutputGermline

  when: 'delly' in tools && runGermline

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idNormal}_${svType}.bcf \
    ${bamNormal}

  delly filter \
    --filter germline \
    --outfile ${idNormal}_${svType}.filter.bcf \
    ${idNormal}_${svType}.bcf

  bcftools view --output-type z ${idNormal}_${svType}.filter.bcf > ${idNormal}_${svType}.delly.vcf.gz
  tabix --preset vcf ${idNormal}_${svType}.delly.vcf.gz
  """
}