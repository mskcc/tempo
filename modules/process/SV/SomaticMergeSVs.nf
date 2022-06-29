process SomaticMergeSVs {
  tag "${idTumor}__${idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.merged.vcf.{gz,gz.tbi}"
  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.merged.raw.vcf"

  input:
    tuple val(idTumor), val(idNormal), val(target), 
      path(Vcfs), path(Tbis),
      val(callerNames)
    path(custom_scripts) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.merged.vcf.gz"), path("${outputPrefix}.merged.vcf.gz.tbi"), emit: SVCallsCombinedVcf
    path("${outputPrefix}.merged.clean.vcf")

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  vcfMap = [:]
  for (i in 1..callerNames.size()){
    vcfMap.put(callerNames[i-1], Vcfs[i-1])
  }
  labelparam = callerNames.sort().join(",")
  inVCFs = ""
  for (i in callerNames.sort()){
    inVCFs += " " + vcfMap[i]
  }
  passMin = callerNames.size() > 2 ? 2 : 1
  """
  mergesvvcf \\
    -n -m ${passMin} \\
    -l ${labelparam} \\
    -o ${outputPrefix}.merged.raw.vcf \\
    -f -d -s -v \\
    ${inVCFs}

  cat ${outputPrefix}.merged.raw.vcf | \\
    awk -F"\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{print}' | \\
  bcftools sort --temp-dir ./ \\
    > ${outputPrefix}.merged.clean.anon.vcf

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
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
