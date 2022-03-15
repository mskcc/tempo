process SomaticAggregateMetadata {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(idTumors), val(idNormals), val(cohort), val(placeHolder), path(metaDataFile)
    
  output:
    path("sample_data.txt"), emit: MetaDataAggregatedOutput

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge metadata file
  mkdir sample_data_tmp
  mv *.sample_data.txt sample_data_tmp/
  awk 'FNR==1 && NR!=1{next;}{print}' sample_data_tmp/*.sample_data.txt > sample_data.txt
  """
}