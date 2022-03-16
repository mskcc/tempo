process GermlineSVVcf2Bedpe {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: ".combined.filtered.bedpe"
  publishDir "${params.outDir}/germline/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.annotsv.tsv"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(vcfFile), path(tbiFile)
    path(repeatMasker)
    path(mapabilityBlacklist)
    path(annotSVref)
    path(spliceSites)
    path(custom_scripts)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.filtered.bedpe"), emit: GermlineSVCombinedBedpe
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.combined.filtered.pass.bedpe"), emit: GermlineSVCombinedBedpePass
    path("${outputPrefix}.combined.bedpe"), emit: GermlineCombinedUnfilteredBedpe
    path("${outputPrefix}.annotsv.tsv"), emit: GermlineCombinedAnnotSV

  when: workflow.profile != "test" 

  script:
  outputPrefix = "${idNormal}"
  genomeBuild = "GRCh37"
  """
  export LC_ALL=C
  /opt/annotsv/bin/AnnotSV \\
    -annotationsDir ${annotSVref} \\
    -SVinputFile ${vcfFile} \\
    -bcftools \$(which bcftools) \\
    -genomeBuild ${genomeBuild} \\
    -includeCI 0 \\
    -outputFile ${outputPrefix}.annotsv

  if [ -f *_AnnotSV/${outputPrefix}.annotsv.tsv ] ; then 
    cp *_AnnotSV/${outputPrefix}.annotsv.tsv .
  else
    touch ${outputPrefix}.annotsv.tsv
  fi
    
  svtools vcftobedpe \\
    -i ${vcfFile} \\
    -o ${outputPrefix}.combined.bedpe \\
    -t ${outputPrefix}_tmp

  if [ ! -s ${outputPrefix}.combined.bedpe ] ; then 
    echo -e "#CHROM_A\\tSTART_A\\tEND_A\\tCHROM_B\\tSTART_B\\tEND_B\\tID\\tQUAL\\tSTRAND_A\\tSTRAND_B\\tTYPE\\tFILTER\\tNAME_A\\tREF_A\\tALT_A\\tNAME_B\\tREF_B\\tALT_B\\tINFO_A\\tINFO_B\\tFORMAT\\t${idNormal}" >> ${outputPrefix}.combined.bedpe
  fi

  python ${custom_scripts}/pair_to_bed_annot.py \\
    --blacklist-regions ${mapabilityBlacklist} \\
    --bedpe ${outputPrefix}.combined.bedpe \\
    --tag mappability \\
    --output ${outputPrefix}.combined.dac.bedpe \\
    --match-type either 
  
  python ${custom_scripts}/pair_to_bed_annot.py \\
    --blacklist-regions ${repeatMasker} \\
    --bedpe ${outputPrefix}.combined.dac.bedpe \\
    --tag repeat_masker \\
    --output ${outputPrefix}.combined.dac.rm.bedpe \\
    --match-type either 
  
  python ${custom_scripts}/detect_cdna.py \\
    --intermediate \\
    --exon-junct ${spliceSites} \\
    --bedpe ${outputPrefix}.combined.dac.rm.bedpe \\
    --out-bedpe ${outputPrefix}.combined.dac.rm.cdna.bedpe \\
    --out ${outputPrefix}.contamination.tsv

  grep -v "^##" ${outputPrefix}.combined.bedpe | \\
    cut -f 1-3,7,11 > ${outputPrefix}.combined.bp1.bed
  echo -n "#" > ${outputPrefix}.combined.bp2.bed
  grep -v "^##" ${outputPrefix}.combined.bedpe | \\
    cut -f 4-6,7,11 >> ${outputPrefix}.combined.bp2.bed
  
  /opt/annotsv/bin/AnnotSV \\
    -annotationsDir ${annotSVref} \\
    -SVinputFile ${outputPrefix}.combined.bp1.bed \\
    -bcftools \$(which bcftools) \\
    -genomeBuild ${genomeBuild} \\
    -includeCI 0 \\
    -svtBEDcol 5 \\
    -outputFile ${outputPrefix}.bp1.annotsv
  /opt/annotsv/bin/AnnotSV \\
    -annotationsDir ${annotSVref} \\
    -SVinputFile ${outputPrefix}.combined.bp2.bed \\
    -bcftools \$(which bcftools) \\
    -genomeBuild ${genomeBuild} \\
    -includeCI 0 \\
    -svtBEDcol 7 \\
    -outputFile ${outputPrefix}.bp2.annotsv

  python ${custom_scripts}/annotate_bedpe.py \\
    --bp1 *_AnnotSV/${outputPrefix}.bp1.annotsv.tsv \\
    --bp2 *_AnnotSV/${outputPrefix}.bp2.annotsv.tsv \\
    --bedpe ${outputPrefix}.combined.dac.rm.cdna.bedpe \\
    --out ${outputPrefix}.combined.dac.rm.cdna.annot.bedpe

  svtools bedpesort \\
    ${outputPrefix}.combined.dac.rm.cdna.annot.bedpe \\
    ${outputPrefix}.combined.filtered.bedpe

  awk -F"\\t" '\$1 ~ /#/ || \$12 == "PASS"' ${outputPrefix}.combined.filtered.bedpe > \\
    ${outputPrefix}.combined.filtered.pass.bedpe
  
  """

}