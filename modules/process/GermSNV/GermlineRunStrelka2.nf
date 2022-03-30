process GermlineRunStrelka2 {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/strelka2", mode: params.publishDirMode

  input:
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal), path(targets), path(targetsIndex)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)
    
  output:
    tuple val("placeHolder"), val(idNormal), val(target), path("${idNormal}.strelka2.vcf.gz"), path("${idNormal}.strelka2.vcf.gz.tbi"), emit: strelkaOutputGermline
  
  script:
  options = ""
  intervals = targets
  if (params.assayType == "exome") {
    options = "--exome"
  }
  """
  configureStrelkaGermlineWorkflow.py \
    ${options} \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Strelka/results/variants/variants.vcf.gz ${idNormal}.strelka2.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi ${idNormal}.strelka2.vcf.gz.tbi
  """
}
