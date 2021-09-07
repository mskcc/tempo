process QcCollectHsMetrics {
  tag {idSample}

  publishDir "${params.outDir}/bams/${idSample}/collecthsmetrics", mode: params.publishDirMode

  input:
    tuple val(idSample), val(target), path(bam), path(bai), path(targetsList), path(baitsList)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)
    val(assayType)
    val(runQC)

  output:
    tuple val(idSample), path("${idSample}.hs_metrics.txt"), emit: collectHsMetricsOutput

  when: params.assayType == "exome" && runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  memMultiplier = params.mem_per_core ? task.cpus : 1
  javaOptions = "--java-options '-Xmx" + task.memory.toString().split(" ")[0].toInteger() * memMultiplier + "g'"

  baitIntervals = "${baitsList}"
  targetIntervals = "${targetsList}"
  """
  gatk CollectHsMetrics \
    ${javaOptions} \
    --TMP_DIR ./ \
    --INPUT ${bam} \
    --OUTPUT ${idSample}.hs_metrics.txt \
    --REFERENCE_SEQUENCE ${genomeFile} \
    --BAIT_INTERVALS ${baitIntervals} \
    --TARGET_INTERVALS ${targetIntervals}
  """
}