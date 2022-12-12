process GermlineSVVcf2Bedpe {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${outputPrefix}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.combined.bedpe"

  input:
    tuple val(idNormal), val(target), path(vcfFile), path(tbiFile)

  output:
    tuple val(idNormal), val(target), path("${outputPrefix}.combined.bedpe"), emit: GermlineCombinedUnfilteredBedpe

  script:
  outputPrefix = "${idNormal}"
  """
  export LC_ALL=C
    
  echo -e "${idNormal} NORMAL" > normalize.samplenames.tsv
  bcftools reheader \\
    --samples normalize.samplenames.tsv \\
    --output reheader_${vcfFile} \\
    ${vcfFile}

  svtools vcftobedpe \\
    -i reheader_${vcfFile} \\
    -o ${outputPrefix}.combined.tmp.bedpe \\
    -t ${outputPrefix}_tmp

  zgrep "^##" ${vcfFile} | sed "s/##fileformat=*/##fileformat=BEDPE/g" > ${outputPrefix}.combined.unsorted.bedpe
  grep -v "^##" ${outputPrefix}.combined.tmp.bedpe | \\
  awk -F"\\t" -v nid="${idNormal}" 'NR == 1 {print \$0,"NORMAL_ID";next;}{print \$0,nid}' \\
  >> ${outputPrefix}.combined.unsorted.bedpe

  if [ ! -s ${outputPrefix}.combined.unsorted.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\tNORMAL\\tNORMAL_ID" >> ${outputPrefix}.combined.unsorted.bedpe
  fi

  svtools bedpesort \\
    ${outputPrefix}.combined.unsorted.bedpe \\
    ${outputPrefix}.combined.bedpe
  """

}
