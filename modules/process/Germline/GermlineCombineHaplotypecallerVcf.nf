process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/haplotypecaller", mode: params.publishDirMode

  input:
    tuple val(id), val(idNormal), val(target), path(haplotypecallerSnpVcf), path(haplotypecallerSnpVcfIndex), path(haplotypecallerIndelVcf), path(haplotypecallerIndelVcfIndex)
    tuple file(genomeFile), file(genomeIndex), file(genomeDict)
    val(tools)
    val(runGermline)

  output:
    tuple val("placeHolder"), val(idNormal), val(target), path("${outfile}"), path("${outfile}.tbi"), emit: haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools && runGermline 

  script: 
  idNormal = id.toString().split("@")[0]
  target = id.toString().split("@")[1]
  outfile = "${idNormal}.haplotypecaller.vcf.gz"  
  """
  bcftools concat \
    --allow-overlaps \
    ${haplotypecallerSnpVcf} ${haplotypecallerIndelVcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all \
    --output-type z \
    --output ${idNormal}.haplotypecaller.vcf.gz

  tabix --preset vcf ${idNormal}.haplotypecaller.vcf.gz
  """
}