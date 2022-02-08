process GermlineRunManta {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/manta", mode: params.publishDirMode

  input:
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex)
    tuple path(svCallingIncludeRegions), path(svCallingIncludeRegionsIndex)

  output:
    tuple val(idNormal), val(target), path("${idNormal}.manta.vcf.gz"), path("${idNormal}.manta.vcf.gz.tbi"), emit: mantaOutputGermline


  // flag with --exome if exome
  script:
  options = ""
  if (params.assayType == "exome") options = "--exome"
  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --reference ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    ${idNormal}.manta.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    ${idNormal}.manta.vcf.gz.tbi
  """
}