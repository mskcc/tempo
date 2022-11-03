process SomaticAggregateHRDetect {
tag { cohortID }
publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

input:
  tuple val(cohortID), 
    path(hrdetectFiles)
output:
  path("hrdetect.tsv")

script:
  """
  awk 'FNR==1 && NR!=1{next;}{print}' ${hrdetectFiles} > hrdetect.tsv
  """
}
