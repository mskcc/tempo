process RunMutationSignatures {
  tag "${idTumor + "__" + idNormal}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(maf)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.mutsig.txt"), emit: mutSig4MetaDataParser

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  maf2cat2.R ${outputPrefix}.somatic.maf \
  ${outputPrefix}.trinucmat.txt
  tempoSig.R --cosmic_${params.cosmic} --pvalue --nperm 10000 --seed 132 ${outputPrefix}.trinucmat.txt \
  ${outputPrefix}.mutsig.txt
  """
}
