process GermlineAggregateSv {
  tag { cohortID }
  publishDir "${params.outDir}/cohort_level/${cohortID}", mode: params.publishDirMode

input:
  tuple val(cohortID),
    path(bedpeFiles)

output:
  path("sv_germline.bedpe") 

script:
  """
  awk '\$1 ~ /^#/ && FNUM !=1{next;}{print}' ${bedpeFiles} > sv_germline.bedpe
  """
}
