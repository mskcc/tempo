process QcPileup {
  tag "${idSample}"

  publishDir "${params.outDir}/bams/${idSample}/pileup/", mode: params.publishDirMode

  input:
    tuple val(idSample), val(target), path(bam), path(bai)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)

  output:
    tuple val(idSample), path("${idSample}.pileup"), emit: pileupOutput

  script:
  gatkPath = "/usr/bin/GenomeAnalysisTK.jar"
  conpairPath = "/usr/bin/conpair"
  markersBed = ""
  if (params.genome == "GRCh37") {
    markersBed = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed"
  }
  else {
    markersBed = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
  }

  if (params.mem_per_core) {
    mem = task.memory.toString().split(" ")[0].toInteger() - 1
  }
  else {
    mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
  }
  javaMem = "${mem}g"
  """
  ${conpairPath}/scripts/run_gatk_pileup_for_sample.py \
    --gatk=${gatkPath} \
    --bam=${bam} \
    --markers=${markersBed} \
    --reference=${genomeFile} \
    --xmx_java=${javaMem} \
    --outfile=${idSample}.pileup
  """
}
