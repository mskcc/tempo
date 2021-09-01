process RunPolysolver {
  tag {idNormal}
  
  input:
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal)
    val(tools)
    val(runSomatic)

  output:
    tuple val("placeHolder"), val(idNormal), val(target), path("${outputPrefix}.hla.txt"), emit: hlaOutput

  when: "polysolver" in tools && runSomatic && ["GRCh38","GRCh37"].contains(params.genome)
  
  script:
  outputPrefix = "${idNormal}"
  outputDir = "."
  tmpDir = "${outputDir}-nf-scratch"
  genome_ = params.genome == "GRCh37" ? "hg19" : "hg38"
  """
  cp /home/polysolver/scripts/shell_call_hla_type .
  
  sed -i "171s/TMP_DIR=.*/TMP_DIR=${tmpDir}/" shell_call_hla_type 

  bash shell_call_hla_type \
    ${bamNormal} \
    Unknown \
    1 \
    ${genome_} \
    STDFQ \
    0 \
    ${outputDir}

  mv winners.hla.txt ${outputPrefix}.hla.txt
  """
}