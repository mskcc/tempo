params.outDir = ""

process GermlineFacetsAnnotation {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/combined_mutations/", mode: params.publishDirMode, pattern: "*.germline.final.maf"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(purity_rdata), path(purity_cncf), path(hisens_cncf), val(facetsPath), path(maf)
    val(tools)
    val(runGermline)
    
  output:
    path("${outputPrefix}.final.maf"), emit: mafFileOutputGermline
    tuple val("placeHolder"), val(idTumor), val(idNormal), file("${outputPrefix}.final.maf"), emit: mafFile4AggregateGermline

  when: tools.containsAll(["facets", "haplotypecaller", "strelka2"]) && runGermline

  script:
  outputPrefix = "${idTumor}__${idNormal}.germline"
  """
  if [ \$( cat ${maf} | wc -l ) -gt 1 ] ; then 
  Rscript --no-init-file /usr/bin/facets-suite/annotate-maf-wrapper.R \
    --facets-output ${purity_rdata} \
    --maf-file ${maf} \
    --output ${outputPrefix}.facets.maf

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-germline.R ${outputPrefix}.facets.maf ${outputPrefix}.final.maf
  else 
    cp ${maf} ${outputPrefix}.final.maf
  fi
  """
}