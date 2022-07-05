process SomaticSVVcf2Bedpe {
  tag "${idTumor}__${idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.combined.bedpe"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfFile), path(tbiFile)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.bedpe"), emit: SomaticCombinedUnfilteredBedpe

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  export LC_ALL=C
    
  svtools vcftobedpe \\
    -i ${vcfFile} \\
    -o ${outputPrefix}.combined.tmp.bedpe \\
    -t ${outputPrefix}_tmp

  zgrep "^#" ${vcfFile} | sed "s/##fileformat=*/##fileformat=BEDPE/g" > ${outputPrefix}.combined.unsorted.bedpe
  grep -v "^#" ${outputPrefix}.combined.tmp.bedpe >> ${outputPrefix}.combined.unsorted.bedpe

  if [ ! -s ${outputPrefix}.combined.unsorted.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\t${idNormal}\\t${idTumor}" >> ${outputPrefix}.combined.unsorted.bedpe
  fi

  svtools bedpesort \\
    ${outputPrefix}.combined.unsorted.bedpe \\
    ${outputPrefix}.combined.bedpe
  """

}
