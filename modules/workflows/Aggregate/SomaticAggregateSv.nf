process SomaticAggregateSv {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(dellyMantaVcf), path(dellyMantaVcfTbi) 
    
  output:
    path("sv_somatic.vcf.{gz,gz.tbi}"), emit: svAggregatedOutput

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge Delly and Manta VCFs
  mkdir sv/
  mv *delly.manta.vcf.gz* sv/
  vcfs=(\$(ls sv/*delly.manta.vcf.gz))
  if [[ \${#vcfs[@]} > 1 ]]
  then
    bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output sv_somatic.vcf.gz \
    sv/*delly.manta.vcf.gz
  else
    mv \${vcfs[0]} sv_somatic.vcf.gz
  fi

  tabix --preset vcf sv_somatic.vcf.gz
  """
}