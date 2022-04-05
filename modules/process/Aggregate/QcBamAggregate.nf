process QcBamAggregate {
  tag "${cohort}"

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(alfredIgnoreYTumor), path(alfredIgnoreYNormal), path(alfredIgnoreNTumor), path(alfredIgnoreNNormal), file(hsMetricsTumor), file(hsMetricsNormal)
    
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
  Rscript --no-init-file /usr/bin/create-aggregate-qc-file.R -n ${task.cpus} -a ${assayType}
  """
}
