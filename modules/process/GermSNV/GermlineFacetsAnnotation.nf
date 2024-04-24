process GermlineFacetsAnnotation {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/combined_mutations/", mode: params.publishDirMode, pattern: "*.germline.final.maf"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(hisens_rdata), val(facetsPath), path(maf)
    
  output:
    path("${outputPrefix}.final.maf"), emit: mafFileOutputGermline
    tuple val(idTumor), val(idNormal), file("${outputPrefix}.final.maf"), emit: mafFile4AggregateGermline

  script:
  outputPrefix = "${idTumor}__${idNormal}.germline"
  """
  if [ \$( cat ${maf} | wc -l ) -gt 1 ] ; then 
  Rscript --no-init-file /usr/bin/facets-suite/annotate-maf-wrapper.R \
    --facets-output ${hisens_rdata} \
    --maf-file ${maf} \
    --output ${outputPrefix}.facets.maf

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-germline.R ${outputPrefix}.facets.maf ${outputPrefix}.final.maf
  else 
    cp ${maf} ${outputPrefix}.final.maf
  fi
  """
}
