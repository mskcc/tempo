process HRDetect {

  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/hrdetect/", mode: params.publishDirMode, pattern: "*.hrdetect.tsv"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(mafFile), path(cnvFile), path(svFile)
    path(HRDetect_script) 

  output:
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("${outputPrefix}.hrdetect.tsv")

  when: params.assayType == "genome" 

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  genome_version = params.genome == 'GRCh38' ? "hg38" : "hg19"
  """
  echo -e "sample\\tsv\\tmutations\\tcnv" > ${outputPrefix}.tsv
  echo -e "${outputPrefix}\\t${svFile}\\t${mafFile}\\t${cnvFile}" >> ${outputPrefix}.tsv
  Rscript ${HRDetect_script} ${outputPrefix}.tsv ${genome_version} ${task.cpus}
  """

}