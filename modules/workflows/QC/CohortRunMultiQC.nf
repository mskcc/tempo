process CohortRunMultiQC {
  tag {cohort}
  label 'multiqc_process'

  publishDir "${params.outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    tuple val(cohort), path(fastPTumor), path(fastPNormal), path(alfredIgnoreYTumor), path(alfredIgnoreYNormal), path(alfredIgnoreNTumor), path(alfredIgnoreNNormal), path(concordFile), path(contamiFile), path(FacetsSummaryFile), path(FacetsQCFile), path(qualimapFolderTumor), path(qualimapFolderNormal), path(hsMetricsTumor), path(hsMetricsNormal)
    tuple path("exome_multiqc_config.yaml"), path("wgs_multiqc_config.yaml"), path("tempoLogo.png")
    
  output:
    tuple val(cohort), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: cohort_multiqc_report

  script:
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
  cp conpair.tsv conpair_genstat.tsv

  mkdir -p fastp_original 
  for i in `find . -maxdepth 1 -name "*fastp.json"` ; do 
    mv \$i fastp_original 
    inname=fastp_original/\$(basename \$i)
    clean_fastp.py \$inname \$i
  done
  
  for i in `find qualimap -name genome_results.txt` ; do
    sampleName=\$(dirname \$i | xargs -n 1 basename )
    cover=\$(grep -i "mean cover" \$i | cut -f 2 -d"=" | sed "s/\\s*//g" | tr -d "X" | tr -d "," )
    echo -e "\${sampleName}\\t\${cover}"
  done > flatCoverage
  echo -e "\\tTumor_Coverage\\tNormal_Coverage" > coverage_split.txt
  join -1 2 -2 1 -o 1.1,1.2,1.3,2.2 -t \$'\\t' <(join -1 1 -2 1 -t \$'\\t' <(cut -f2,3 conpair_genstat.tsv | tail -n +2 | sort | uniq) <(cat flatCoverage | sort | uniq)) <(cat flatCoverage | sort | uniq) | cut -f 1,3,4 >> coverage_split.txt
  join -1 1 -2 1 -t \$'\\t' <(cut -f 3 conpair_genstat.tsv | sort | uniq | sed "s/\$/\\t/g" ) <(cat flatCoverage | sort | uniq) >> coverage_split.txt
  
  parse_alfred.py --alfredfiles *alfred*tsv.gz
  mkdir -p ignoreFolder 
  find . -maxdepth 1 \\( -name 'CO_ignore_mqc.yaml' -o -name 'IS_*mqc.yaml' -o -name 'GC_ignore_mqc.yaml' -o -name 'ME_aware_mqc.yaml' \\) -type f -print0 | xargs -0r mv -t ignoreFolder
  if [[ "${params.assayType}" == "exome" ]] ; then 
    find . -maxdepth 1 -name 'CM_*mqc.yaml' -type f -print0 | xargs -0r mv -t ignoreFolder
  fi
  mv conpair.tsv ignoreFolder

  for i in *.facets_qc.txt ; do 
    head -1 \$i | cut -f 1,28,97  | sed "s/^tumor_sample_id//g"> \$i.qc.txt
    tail -n +2 \$i | cut -f 1,28,97 | sed "s/TRUE\$/PASS/g" | sed "s/FALSE\$/FAIL/g" >> \$i.qc.txt
  done

  cp ${assay}_multiqc_config.yaml multiqc_config.yaml

  samplesNum=`for i in ./*contamination.txt ; do tail -n +2 \$i | cut -f 2 ; done | sort | uniq | wc -l`
  fastpNum=`ls ./*fastp*json | wc -l`
  mqcSampleNum=\$((samplesNum + fastpNum ))
  mqcSampleNum=\$(( mqcSampleNum > 25 ? mqcSampleNum : 25 ))
  
  multiqc . --cl_config "max_table_rows: \$(( mqcSampleNum + 1 ))" -x ignoreFolder/ -x fastp_original/
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  if [ \$samplesNum -gt 50 ] ; then 
    cp genstats-QC_Status.txt QC_Status.txt
  fi
  beeswarm_config="max_table_rows: \$(( mqcSampleNum + 1 ))"

  multiqc . --cl_config "title: \\"Cohort MultiQC Report\\"" --cl_config "subtitle: \\"${cohort} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "\${beeswarm_config}" --cl_config "report_comment: \\"This report includes FASTQ and alignment for all samples in ${cohort}and Tumor/Normal pair statistics for all pairs in ${cohort}.\\"" -z -x ignoreFolder/ -x fastp_original/

  """
}
