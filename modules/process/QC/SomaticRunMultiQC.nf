process SomaticRunMultiQC {
   tag {idTumor + "__" + idNormal}
   label 'multiqc_process'

  publishDir "${params.outDir}/somatic/${outPrefix}/multiqc", mode: params.publishDirMode  
  
  input:
    tuple val(idTumor), val(idNormal), path(conpairFiles), path(facetsSummaryFiles), path(facetsQCFiles)
    tuple path("exome_multiqc_config.yaml"), path("wgs_multiqc_config.yaml"), path("tempoLogo.png")

  output:
    tuple val(idTumor), val(idNormal), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: somatic_multiqc_report

  script: 
  outPrefix = "${idTumor}__${idNormal}"
  if (params.assayType == "exome") {
    assay = "exome"
  }
  else {
    assay = 'wgs'
  }
  """
  echo -e "\\tTumor\\tNormal\\tTumor_Contamination\\tNormal_Contamination\\tConcordance" > conpair.tsv
  for i in ./*contamination.txt ; do 
     j=./\$(basename \$i | cut -f 1 -d.).concordance.txt
     echo -e "\$(tail -n +2 \$i | sort -r | cut -f 2| head -1)\\t\$(tail -n +2 \$i | sort -r | cut -f 2| paste -sd"\\t")\\t\$(tail -n +2 \$i | sort -r | cut -f 3| paste -sd"\\t")\\t\$(tail -1 \$j | cut -f 2 )" >> conpair.tsv
  done

  head -1 ${facetsQCFiles} | cut -f 1,28,97  | sed "s/^tumor_sample_id//g"> ${facetsQCFiles}.qc.txt
  tail -n +2 ${facetsQCFiles} | cut -f 1,28,97 | sed "s/TRUE\$/PASS/g" | sed "s/FALSE\$/FAIL/g" >> ${facetsQCFiles}.qc.txt

  cp conpair.tsv conpair_genstat.tsv
  mkdir -p ignoreFolder ; mv conpair.tsv ignoreFolder
  cp ${assay}_multiqc_config.yaml multiqc_config.yaml
  
  multiqc . -x ignoreFolder
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  multiqc . --cl_config "title: \\"Somatic MultiQC Report\\"" --cl_config "subtitle: \\"${outPrefix} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "report_comment: \\"This report includes QC statistics related to the Tumor/Normal pair ${outPrefix}.<br/>This report does not include FASTQ or alignment QC of either ${idTumor} or ${idNormal}. To review FASTQ and alignment QC, please refer to the multiqc_report.html from the bam-level folder.<br/>To review qc from all samples and Tumor/Normal pairs from a cohort in a single report, please refer to the multiqc_report.html from the cohort-level folder.\\"" -z -x ignoreFolder

  """

}