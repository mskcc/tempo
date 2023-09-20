process RunPlatypus {
  tag "${idTumor + "__" + idNormal}"
  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/platypus", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict) 

  output:
    tuple val(combineKey), path('*Somatic*'), emit:  platypusCombine
    path("${outputPrefix}.Somatic.Platypus.vcf"), emit: platypusOutput

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  outfile = "${outputPrefix}.Somatic.Platypus.vcf"
  combineKey = "${idTumor}__${idNormal}_${target}"

  """

    python /code/Platypus/bin/Platypus.py callVariants \
    --bamFiles=${bamNormal},${bamTumor} \
    --refFile=${genomeFile} \
    --output=${outputPrefix}_bothPlat.vcf \
    --nCPU=2 \
    --genSNPs=1 \
    --minVarFreq=.05

    python /code/Platypus/extensions/Cancer/somaticMutationDetector.py \
    --inputVCF ${outputPrefix}_bothPlat.vcf \
    --outputVCF ${outputPrefix}.Somatic.Platypus.vcf \
    --tumourSample  ${idTumor} \
    --normalSample ${idNormal}

  """
}
