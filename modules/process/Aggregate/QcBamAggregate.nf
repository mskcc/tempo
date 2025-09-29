process QcBamAggregate {
  tag "${cohort}"

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(alfredIgnoreYTumor), path(alfredIgnoreYNormal), path(alfredIgnoreNTumor), path(alfredIgnoreNNormal), file(hsMetricsTumor), file(hsMetricsNormal)
    path(aggregate_qc_Rscript)

  output:
    path('alignment_qc.txt'), emit: alignmentQcAggregatedOutput

  script:
  if (params.assayType == "exome") {
    assayType = "wes"
  }
  else {
    assayType = 'wgs'
  }
  """
  Rscript --no-init-file ${aggregate_qc_Rscript} -n ${task.cpus} -a ${assayType}
  """
}
