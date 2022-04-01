process GermlineSVVcf2Bedpe {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.combined.bedpe"

  input:
    tuple val(idNormal), val(target), path(vcfFile), path(tbiFile)

  output:
    tuple val(idNormal), val(target), path("${outputPrefix}.combined.bedpe"), emit: GermlineCombinedUnfilteredBedpe

  script:
  outputPrefix = "${idNormal}"
  """
  export LC_ALL=C
    
  svtools vcftobedpe \\
    -i ${vcfFile} \\
    -o ${outputPrefix}.combined.unsorted.bedpe \\
    -t ${outputPrefix}_tmp

  if [ ! -s ${outputPrefix}.combined.unsorted.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\t${idNormal}" >> ${outputPrefix}.combined.unsorted.bedpe
  fi

  svtools bedpesort \\
    ${outputPrefix}.combined.unsorted.bedpe \\
    ${outputPrefix}.combined.bedpe
  
  """

}