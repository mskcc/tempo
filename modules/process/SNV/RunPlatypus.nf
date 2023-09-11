process RunPlatypus {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/platypus", mode: params.publishDirMode, pattern: "*.vcf.{gz,gz.tbi}"

  input:
    tuple val(id), val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(intervalBed)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path('*Platypus.vcf'), emit:  platypusCombine
    path('*Somatic.vcf'), emit: platypusOutput

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  outfile = "${outputPrefix}.Somatic.Platypus.vcf"
  """

    python /usr/bin/Platypus.py callVariants \
    --bamFiles=${bamNormal},${bamTumor} \
    --refFile=${genomeFile} \
    --output=${outputPrefix}_bothPlat.vcf \
    --nCPU=2 \
    --genSNPs=1 \
    --minVarFreq=.05

    python /usr/extensions/Cancer/somaticMutationDetector.py \
    --inputVCF ${outputPrefix} \
    --outputVCF ${outfile} \
    --tumourSample  ${idTumor} \
    --normalSample ${idNormal}

  """
}
