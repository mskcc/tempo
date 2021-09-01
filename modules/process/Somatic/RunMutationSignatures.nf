process RunMutationSignatures {
  tag {idTumor + "__" + idNormal}

  input:
    tuple val(idTumor), val(idNormal), val(target), path(maf)
    val(tools)
    val(runSomatic)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.mutsig.txt"), emit: mutSig4MetaDataParser

  when: tools.containsAll(["mutect2", "manta", "strelka2", "mutsig"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  maf2cat2.R ${outputPrefix}.somatic.maf \
  ${outputPrefix}.trinucmat.txt
  tempoSig.R --cosmic_${params.cosmic} --pvalue --nperm 10000 --seed 132 ${outputPrefix}.trinucmat.txt \
  ${outputPrefix}.mutsig.txt
  """
}