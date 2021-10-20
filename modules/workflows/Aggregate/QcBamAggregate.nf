process QcBamAggregate {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(alfredIgnoreYTumor), path(alfredIgnoreYNormal), path(alfredIgnoreNTumor), path(alfredIgnoreNNormal), path(hsMetricsTumor), path(hsMetricsNormal)
    val(runQC)
    
  output:
    path('alignment_qc.txt'), emit: alignmentQcAggregatedOutput

  when: runQC

  script:
  if (params.assayType == "exome") {
    assayType = "wes"
  }
  else {
    assayType = 'wgs'
  }
  """
  Rscript --no-init-file /usr/bin/create-aggregate-qc-file.R -n ${task.cpus} -a ${assayType}
  """
}