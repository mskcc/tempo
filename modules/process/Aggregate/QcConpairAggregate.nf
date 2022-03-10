process QcConpairAggregate {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(concordFile), path(contamiFile)
    
  output:
    tuple path('concordance_qc.txt'), path('contamination_qc.txt'), emit: conpairAggregatedOutput

  script:
  """
  if ls *.concordance.txt 1> /dev/null 2>&1; then
    echo -e "Pair\tConcordance" > concordance_qc.txt
    grep -v "concordance" *.concordance.txt | sed 's/.concordance.txt:/\t/' | cut -f1,3 | sort -k1,1 >> concordance_qc.txt
  fi
  if ls *.contamination.txt 1> /dev/null 2>&1; then
    echo -e "Pair\tSample_Type\tSample_ID\tContamination" > contamination_qc.txt
    grep -v "Contamination" *.contamination.txt | sed 's/.contamination.txt:/\t/' | sort -k1,1 >> contamination_qc.txt
  fi
  touch concordance_qc.txt contamination_qc.txt
  """
}