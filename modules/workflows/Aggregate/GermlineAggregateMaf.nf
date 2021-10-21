process GermlineAggregateMaf {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(idTumors), val(idNormals), val(cohort), val(placeHolder), path(mafFile)

  output:
    path("mut_germline.maf"), emit: mutationAggregatedGermlineOutput

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
  done | grep ^Hugo | head -n1 > mut_germline.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo" | sort -k5,5V -k6,6n >> mut_germline.maf

  """
}
