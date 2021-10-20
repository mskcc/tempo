process GermlineAggregateSv {
  tag {cohort}

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(dellyMantaVcf), path(dellyMantaVcfTbi)
    val(runGermline)
    
  output:
    path("sv_germline.vcf.{gz,gz.tbi}"), emit: svAggregatedGermlineOutput

  when: runGermline

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge Delly and Manta VCFs
  mkdir sv
  mv  *.delly.manta.vcf.gz* sv/
  vcfs=(\$(ls sv/*delly.manta.vcf.gz))
  if [[ \${#vcfs[@]} > 1 ]]
  then
    bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output sv_germline.vcf.gz \
    sv/*delly.manta.vcf.gz
  else
    mv \${vcfs[0]} sv_germline.vcf.gz
  fi

  tabix --preset vcf sv_germline.vcf.gz
  """
}