process SomaticMergeSVs {
  tag "${idTumor}__${idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.merged.vcf.{gz,gz.tbi}"
  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.merged.clean.vcf"

  input:
    tuple val(idTumor), val(idNormal), val(target), 
      path(dellyVcfs), path(dellyVcfsIndex), 
      path(mantaFile), path(mantaIndex), 
      file(svabaFile), file(svabaIndex), 
      file(brassFile), file(brassIndex)
    path(custom_scripts) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.merged.vcf.gz"), path("${outputPrefix}.merged.vcf.gz.tbi"), emit: SVCallsCombinedVcf
    path("${outputPrefix}.merged.clean.vcf")

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  labelparam = "delly,manta" 
  labelparam = svabaFile.name.endsWith("vcf.gz") ?  labelparam + ",svaba" : labelparam 
  labelparam = brassFile.name.endsWith("vcf.gz") ?  labelparam + ",brass" : labelparam 
  labelparam_list = labelparam.split(",")
  inVCFs = "${outputPrefix}.delly.vcf.gz ${mantaFile} "
  inVCFs = labelparam_list.contains("svaba") ? inVCFs + " svaba.reformat.vcf.gz " : inVCFs 
  inVCFs = labelparam_list.contains("brass") ? inVCFs + " brass.reformat.vcf.gz " : inVCFs 
  passMin = labelparam_list.size() > 2 ? 2 : 1
  """
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${outputPrefix}.delly.vcf.gz \\
    *.delly.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.vcf.gz 

  if [ "${labelparam_list.contains("brass")}" = true ] ; then 
    echo -e "TUMOUR ${idTumor}\\nNORMAL ${idNormal}" > brass.samplenames.tsv 
    bcftools reheader \\
      --samples brass.samplenames.tsv \\
      --output brass.reformat.vcf.gz \\
      ${brassFile}
    tabix -p vcf brass.reformat.vcf.gz
  fi

  if [ "${labelparam_list.contains("svaba")}" = true ] ; then 
    echo -e "${idTumor}.bam ${idTumor}\\n${idNormal}.bam ${idNormal}" > svaba.samplenames.tsv 
    bcftools reheader \\
      --samples svaba.samplenames.tsv \\
      --output svaba.reformat.vcf.gz \\
      ${svabaFile}
    tabix -p vcf svaba.reformat.vcf.gz
  fi

  mergesvvcf \\
    -n -m ${passMin} \\
    -l ${labelparam} \\
    -o ${outputPrefix}.merged.raw.vcf \\
    -f -d -s -v \\
    ${inVCFs}

  cat ${outputPrefix}.merged.raw.vcf | \\
    awk -F"\\t" -v OFS="\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = "."; print \$0}' | \\
  bcftools sort --temp-dir ./ \\
    > ${outputPrefix}.merged.clean.anon.vcf

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS\\_%INFO/CHR2\\_%INFO/END\\_%INFO/STRANDS' \\
    -o ${outputPrefix}.merged.clean.vcf \\
    ${outputPrefix}.merged.clean.anon.vcf

  python ${custom_scripts}/filter-sv-vcf.py \\
    --input ${outputPrefix}.merged.clean.vcf \\
    --output ${outputPrefix}.merged.clean.corrected.vcf \\
    --min ${passMin}

  bcftools view \\
    --samples ${idTumor},${idNormal} \\
    --output-type z \\
    --output-file ${outputPrefix}.merged.vcf.gz \\
    ${outputPrefix}.merged.clean.corrected.vcf
  
  tabix --preset vcf ${outputPrefix}.merged.vcf.gz 
  """
}
