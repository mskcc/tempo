process RunPolysolver {
  tag {idNormal}
  
  input:
    tuple val(idNormal), val(target), path(bamNormal), path(baiNormal)

  output:
    tuple val("placeHolder"), val(idNormal), val(target), path("${outputPrefix}.hla.txt"), emit: hlaOutput
  
  script:
  outputPrefix = "${idNormal}"
  outputDir = "."
  tmpDir = "${outputDir}-nf-scratch"
  genome_ = params.genome == "GRCh37" ? "hg19" : params.genome == 'GRCh38' ? "hg38" : params.genome == 'smallGRCh37' ? "small" : "other"
  """
if [ ${genome_} != "small" ] ; then

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

else

  echo -e "HLA-A\thla_a_01_01_01_01\thla_a_01_01_01_01\nHLA-B\thla_b_15_02_01\thla_b_15_02_01\nHLA-C\thla_c_01_02_01\thla_c_01_02_01" > ${outputPrefix}.hla.txt

fi
  """
}
