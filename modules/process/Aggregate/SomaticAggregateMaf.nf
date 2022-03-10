process SomaticAggregateMaf {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(idTumors), val(idNormals), val(cohort), val(placeHolder), path(mafFile)

  output:
    path("mut_somatic.maf"), emit: mutationAggregatedOutput

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge MAF files
  mkdir mut
  mv *.maf mut/
  for i in mut/*.maf ; do 
    if [ \$( cat \$i | wc -l ) -gt 1 ] ; then 
      cat \$i
    fi
  done | grep ^Hugo_Symbol | head -n 1 > mut_somatic.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo_Symbol" | sort -k5,5V -k6,6n >> mut_somatic.maf
  """
}