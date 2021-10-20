process RunLOHHLA {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/lohhla", mode: params.publishDirMode

  input:
    tuple val(idNormal), val(target), val(idTumor), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(purityOut), val(placeHolder), path(winnersHla)
    tuple path(hlaFasta), path(hlaDat)

  output:
    tuple path("*.DNA.HLAlossPrediction_CI.txt"), path("*DNA.IntegerCPN_CI.txt"), path("*.pdf"), path("*.RData"), optional: true, emit: lohhlaOutput
    tuple val("placeHolder"), val(idTumor), val(idNormal), file("*.DNA.HLAlossPrediction_CI.txt"), emit: predictHLA4Aggregate
    tuple val("placeHolder"), val(idTumor), val(idNormal), file("*DNA.IntegerCPN_CI.txt"), emit: intCPN4Aggregate

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  cat ${winnersHla} | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
  
  PURITY=\$(grep Purity *_purity.out | grep -oP "[0-9\\.]+|NA+")
  PLOIDY=\$(grep Ploidy *_purity.out | grep -oP "[0-9\\.]+|NA+")
  cat <(echo -e "tumorPurity\ttumorPloidy") <(echo -e "${idTumor}\t\$PURITY\t\$PLOIDY") > tumor_purity_ploidy.txt

  Rscript --no-init-file /lohhla/LOHHLAscript.R \
    --patientId ${outputPrefix} \
    --normalBAMfile ${bamNormal} \
    --tumorBAMfile ${bamTumor} \
    --HLAfastaLoc ${hlaFasta} \
    --HLAexonLoc ${hlaDat} \
    --CopyNumLoc tumor_purity_ploidy.txt \
    --minCoverageFilter ${params.lohhla.minCoverageFilter} \
    --hlaPath massaged.winners.hla.txt \
    --gatkDir /picard-tools \
    --novoDir /opt/conda/bin

  if [[ -f ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt ]]
  then
    sed -i "s/^${idTumor}/${outputPrefix}/g" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt
  else
    rm -rf *.DNA.HLAlossPrediction_CI.txt
  fi

  if [[ -f ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt ]]
  then
    sed -i "s/^/${outputPrefix}\t/g" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt
    sed -i "0,/^${outputPrefix}\t/s//sample\t/" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt
  else
    rm -rf *.DNA.IntegerCPN_CI.txt
  fi

  touch ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt
  touch ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt

  mv ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt ${outputPrefix}.DNA.HLAlossPrediction_CI.txt
  mv ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt ${outputPrefix}.DNA.IntegerCPN_CI.txt

  if find Figures -mindepth 1 | read
  then
    mv Figures/* .
    mv ${idTumor}.minCoverage_${params.lohhla.minCoverageFilter}.HLA.pdf ${outputPrefix}.HLA.pdf
  fi
  """
}