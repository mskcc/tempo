process SomaticSVVcf2Bedpe {
  tag "${idTumor}__${idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.combined.bedpe"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfFile), path(tbiFile)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.bedpe"), emit: SomaticCombinedUnfilteredBedpe

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  export LC_ALL=C

  echo -e "${idTumor} TUMOR\\n${idNormal} NORMAL" > normalize.samplenames.tsv
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
  awk -F"\\t" -v tid="${idTumor}" -v nid="${idNormal}" -v OFS="\\t" 'NR == 1 {print \$0,"TUMOR_ID","NORMAL_ID";next;}{print \$0,tid,nid}' \\
  >> ${outputPrefix}.combined.unsorted.bedpe

  if [ ! -s ${outputPrefix}.combined.unsorted.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\tTUMOR\\tNORMAL\\tTUMOR_ID\\tNORMAL_ID" >> ${outputPrefix}.combined.unsorted.bedpe
  fi

  svtools bedpesort \\
    ${outputPrefix}.combined.unsorted.bedpe \\
    ${outputPrefix}.combined.bedpe
  """

}
