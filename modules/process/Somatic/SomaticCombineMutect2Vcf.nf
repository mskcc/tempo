process SomaticCombineMutect2Vcf {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/mutect2", mode: params.publishDirMode

  input:
    tuple val(id), val(idTumor), val(idNormal), val(target), path(mutect2Vcf), path(mutect2VcfIndex), path(mutect2Stats)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)
    val(tools)
    val(runSomatic)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outfile}"), path("${outfile}.tbi"), emit: mutect2CombinedVcfOutput

  when: "mutect2" in tools && runSomatic

  script:
  idTumor = id.toString().split("__")[0]
  idNormal = id.toString().split("@")[0].split("__")[1]
  target = id.toString().split("@")[1]
  outfile = "${idTumor}__${idNormal}.mutect2.vcf.gz"
  """
  bcftools concat \
    --allow-overlaps \
    ${mutect2Vcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all | \
  bcftools view \
    --samples ${idNormal},${idTumor} \
    --output-type z \
    --output-file ${outfile}

  tabix --preset vcf ${outfile}
  """
}
