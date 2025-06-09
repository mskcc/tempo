process SomaticRunMultiQC {
   tag "${idTumor + "__" + idNormal}"
   label 'multiqc_process'

  publishDir "${params.outDir}/somatic/${outPrefix}/multiqc", mode: params.publishDirMode  
  
  input:
    tuple val(idTumor), val(idNormal), file(conpairFiles), file(qualimapTumor), file(qualimapNormal), file(facetsSummaryFiles), file(facetsQCFiles)
    tuple file("exome_multiqc_config.yaml"), file("wgs_multiqc_config.yaml"), file("tempoLogo.png")

  output:
    tuple val(idTumor), val(idNormal), file("*multiqc_report*.html"), file("*multiqc_data*.zip"), emit: somatic_multiqc_report
    tuple val(idTumor), val(idNormal), file("${outPrefix}.QC_Status.txt")

  script: 
  outPrefix = "${idTumor}__${idNormal}"
  if (params.assayType == "exome") {
    assay = "exome"
  }
  else {
    assay = 'wgs'
  }
  """
  for i in ./*_qualimap_rawdata.tar.gz ; do 
    newFolder=\$(basename \$i | rev | cut -f 3- -d. | cut -f 3- -d_ | rev ) 
    mkdir -p qualimap/\$newFolder
    tar -xzf \$i -C qualimap/\$newFolder
  done

  echo -e "\\tTumor\\tNormal\\tTumor_Contamination\\tNormal_Contamination\\tConcordance" > conpair.tsv
  for i in ./*contamination.txt ; do 
     j=./\$(basename \$i | cut -f 1 -d.).concordance.txt
     echo -e "\$(tail -n +2 \$i | sort -r | cut -f 2| head -1)\\t\$(tail -n +2 \$i | sort -r | cut -f 2| paste -sd"\\t")\\t\$(tail -n +2 \$i | sort -r | cut -f 3| paste -sd"\\t")\\t\$(tail -1 \$j | cut -f 2 )" >> conpair.tsv
  done

  head -1 ${facetsQCFiles} | cut -f 1,28,97  | sed "s/^tumor_sample_id//g"> ${facetsQCFiles}.qc.txt
  tail -n +2 ${facetsQCFiles} | cut -f 1,28,97 | sed "s/TRUE\$/PASS/g" | sed "s/FALSE\$/FAIL/g" >> ${facetsQCFiles}.qc.txt

  cp conpair.tsv conpair_genstat.tsv

  for i in `find qualimap -name genome_results.txt` ; do
    sampleName=\$(dirname \$i | xargs -n 1 basename )
    cover=\$(grep -i "mean cover" \$i | cut -f 2 -d"=" | sed "s/\\s*//g" | tr -d "X" | tr -d "," )
    echo -e "\${sampleName}\\t\${cover}"
  done > flatCoverage
  echo -e "\\tTumor_Coverage\\tNormal_Coverage" > coverage_split.txt
  join -1 2 -2 1 -o 1.1,1.2,1.3,2.2 -t \$'\\t' <(join -1 1 -2 1 -t \$'\\t' <(cut -f2,3 conpair_genstat.tsv | tail -n +2 | sort | uniq) <(cat flatCoverage | sort | uniq)) <(cat flatCoverage | sort | uniq) | cut -f 1,3,4 >> coverage_split.txt
  join -1 1 -2 1 -t \$'\\t' <(cut -f 3 conpair_genstat.tsv | sort | uniq | sed "s/\$/\\t/g" ) <(cat flatCoverage | sort | uniq) >> coverage_split.txt


  mkdir -p ignoreFolder ; cp conpair.tsv ignoreFolder
  cp ${assay}_multiqc_config.yaml multiqc_config.yaml
  
  echo -e "metric\tpass\twarn\tfail\ndummy\tdummy\tdummy\tdummy" > CriteriaTable.txt
  multiqc . -x ignoreFolder -x qualimap
  rm -f CriteriaTable.txt
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  multiqc . --cl_config "title: \\"Somatic MultiQC Report\\"" --cl_config "subtitle: \\"${outPrefix} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "report_comment: \\"This report includes QC statistics related to the Tumor/Normal pair ${outPrefix}.<br/>This report does not include FASTQ or alignment QC of either ${idTumor} or ${idNormal}. To review FASTQ and alignment QC, please refer to the multiqc_report.html from the bam-level folder.<br/>To review qc from all samples and Tumor/Normal pairs from a cohort in a single report, please refer to the multiqc_report.html from the cohort-level folder.\\"" -z -x ignoreFolder -x qualimap
  mv genstats-QC_Status.txt ${outPrefix}.QC_Status.txt

  """

}
