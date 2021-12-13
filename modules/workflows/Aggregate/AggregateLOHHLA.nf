process AggregateLOHHLA {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(preditHLA), path(intCPN)

  output:
    path("DNA.IntegerCPN_CI.txt"), emit: lohhlaDNAIntegerCPNOutput
    path("HLAlossPrediction_CI.txt"), emit: lohhlaHLAlossPredictionAggregatedOutput

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp
  mkdir lohhla
  mv *.txt lohhla/
  awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*HLAlossPrediction_CI.txt > HLAlossPrediction_CI.txt
  awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*DNA.IntegerCPN_CI.txt > DNA.IntegerCPN_CI.txt
  """
}