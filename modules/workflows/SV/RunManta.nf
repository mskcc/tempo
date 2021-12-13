process RunManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/manta", mode: params.publishDirMode, pattern: "*.manta.vcf.{gz,gz.tbi}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex)
    tuple path(svCallingIncludeRegions), path(svCallingIncludeRegionsIndex)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.manta.vcf.gz"), emit: manta4Combine
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.manta.vcf.gz.tbi"), emit: mantaOutput
    tuple val(idTumor), val(idNormal), val(target), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSmallIndels.vcf.gz.tbi"), emit: mantaToStrelka

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  options = ""
  if (params.assayType == "exome") options = "--exome"
  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --referenceFasta ${genomeFile} \
    --normalBam ${bamNormal} \
    --tumorBam ${bamTumor} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${outputPrefix}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${outputPrefix}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${outputPrefix}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    ${outputPrefix}.manta.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    ${outputPrefix}.manta.vcf.gz.tbi
  """
}