process SampleRunMultiQC {
  tag {idSample}
  label 'multiqc_process'

  publishDir "${params.outDir}/bams/${idSample}/multiqc", mode: params.publishDirMode  
  
  input:
    tuple val(idSample), path(alfredRGNTsvFile), path(alfredRGYTsvFile), path(fastpJsonFile), path(qualimapFolder), path(hsmetricsFile)
    tuple path("exome_multiqc_config.yaml"), path("wgs_multiqc_config.yaml"), path("tempoLogo.png")

  output:
    tuple val(idSample), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: sample_multiqc_report
    tuple val(idSample), path("QC_Status.txt")

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
  
  parse_alfred.py --alfredfiles *alfred*tsv.gz 
  mkdir -p ignoreFolder 
  find . -maxdepth 1 \\( -name 'CO_ignore*mqc.yaml' -o -name 'IS_*mqc.yaml' -o -name 'GC_ignore*mqc.yaml' -o -name 'ME_aware_mqc.yaml' \\) -type f -print0 | xargs -0r mv -t ignoreFolder
  if [[ "${params.assayType}" == "exome" ]] ; then 
    find . -maxdepth 1 -name 'CM_*mqc.yaml' -type f -print0 | xargs -0r mv -t ignoreFolder
  fi

  mkdir -p fastp_original 
  for i in `find . -maxdepth 1 -name "*fastp.json"` ; do 
    mv \$i fastp_original 
  done
  for i in `find fastp_original -name "*fastp.json"` ; do
    outname=\$(basename \$i)
    python -c "import json
  with open('\$i', 'r') as data_file:
    data = json.load(data_file)
  data_file.close()
  keys = list(data.keys())
  for element in keys:
    if element in ['read1_after_filtering','read2_after_filtering'] :
      data.pop(element, None)
  keys = list(data['summary'].keys())
  for element in keys:
    if element in ['after_filtering'] :
      data['summary'].pop(element, None)
  with open('\$outname', 'w') as data_file:
    json.dump(data, data_file)
  data_file.close()"
  done

  echo -e "\\tCoverage" > coverage_split.txt
  cover=\$(grep -i "mean cover" ./qualimap/${idSample}/genome_results.txt | cut -f 2 -d"=" | sed "s/\\s*//g" | tr -d "X")
  echo -e "${idSample}\\t\${cover}" >> coverage_split.txt
  
  cp ${assay}_multiqc_config.yaml multiqc_config.yaml

  multiqc . -x ignoreFolder/ -x fastp_original/
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  multiqc . --cl_config "title: \\"Sample MultiQC Report\\"" --cl_config "subtitle: \\"${idSample} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "report_comment: \\"This report includes FASTQ and alignment statistics for the sample ${idSample}.<br/>This report does not include QC metrics from the Tumor/Normal pair that includes ${idSample}. To review pairing QC, please refer to the multiqc_report.html from the somatic-level folder.<br/>To review qc from all samples and Tumor/Normal pairs from a cohort in a single report, please refer to the multiqc_report.html from the cohort-level folder.\\"" -t "tempo" -z -x ignoreFolder/ -x fastp_original/
  mv genstats-QC_Status.txt QC_Status.txt
  """

}