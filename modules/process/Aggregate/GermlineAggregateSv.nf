process GermlineAggregateSv {
  tag "${cohort}"

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(combinedVcf), path(combinedVcfTbi)
    
  output:
    path("sv_germline.vcf.{gz,gz.tbi}"), emit: svAggregatedGermlineOutput

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge Delly and Manta VCFs
  mkdir sv
  mv ${combinedVcf.join(" ")} sv/
  mv ${combinedVcfTbi.join(" ")} sv/
  vcfs=(\$(ls sv/*.vcf.gz))
  if [[ \${#vcfs[@]} > 1 ]]
  then
    bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output sv_germline.vcf.gz \
    sv/*.vcf.gz
  else
    mv \${vcfs[0]} sv_germline.vcf.gz
  fi

  tabix --preset vcf sv_germline.vcf.gz
  """
}
