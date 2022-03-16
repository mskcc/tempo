process SomaticSVVcf2Bedpe {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: ".combined.filtered.bedpe"
  // publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.annotsv.tsv"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfFile), path(tbiFile)

  output:
    //tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.filtered.bedpe"), emit: SVCombinedBedpe
    //tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.filtered.pass.bedpe"), emit: SVCombinedBedpePass
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.bedpe"), emit: SomaticCombinedUnfilteredBedpe
    // path("${outputPrefix}.annotsv.tsv"), emit: SomaticCombinedAnnotSV

  when: workflow.profile != "test" 

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  genomeBuild = "GRCh37"
  """
  export LC_ALL=C
    
  svtools vcftobedpe \\
    -i ${vcfFile} \\
    -o ${outputPrefix}.combined.unsorted.bedpe \\
    -t ${outputPrefix}_tmp

  if [ ! -s ${outputPrefix}.combined.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\t${idNormal}\\t${idTumor}" >> ${outputPrefix}.combined.unsorted.bedpe
  fi

  svtools bedpesort \\
    ${outputPrefix}.combined.unsorted.bedpe \\
    ${outputPrefix}.combined.bedpe
  """

}
