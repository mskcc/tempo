params.outDir = ""

process DoFacetsPreviewQC {
  tag {idTumor + "__" + idNormal}
  publishDir "${params.outDir}/somatic/${tag}/facets/${tag}/", mode: params.publishDirMode, pattern: "${idTumor}__${idNormal}.facets_qc.txt"

  input:
    tuple val(idTumor), val(idNormal), val(target), file(facetsOutputFolderFiles), path(countsFile), val(facetsOutputDir)
    val(tools)
    val(runSomatic)

  output:
    tuple val(idTumor), val(idNormal), path("${idTumor}__${idNormal}.facets_qc.txt"), emit: FacetsPreviewOut

  when: "facets" in tools && runSomatic

  script:
  tag = "${idTumor}__${idNormal}"
  """
  mkdir -p ${facetsOutputDir} 
  facetsFitFiles=( ${facetsOutputFolderFiles.join(" ")} )
  for i in "\${facetsFitFiles[@]}" ; do 
    cp \$i ${facetsOutputDir}/\$i
  done
  echo -e "sample_id\\tsample_path\\ttumor_id" > manifest.txt 
  echo -e "${idTumor}__${idNormal}\\t\$(pwd)\\t${idTumor}" >> manifest.txt 
  gzip manifest.txt
  mkdir -p refit_watcher/bin/ refit_watcher/refit_jobs/
  R -e "facetsPreview::generate_genomic_annotations('${idTumor}__${idNormal}', '\$(pwd)/', '/usr/bin/facets-preview/tempo_config.json')"
  cp facets_qc.txt ${idTumor}__${idNormal}.facets_qc.txt
  rm ${facetsOutputDir}/*
  """

}