process GermlineMergeSVs {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/", mode: params.publishDirMode, pattern: "*.merged.vcf.{gz,gz.tbi}"

  input:
    tuple val(idNormal), val(target), 
      path(Vcfs), path(Tbis),
      path(callerNames)
    path(custom_scripts)

  output:
    tuple val(idNormal), val(target), path("${idNormal}.merged.vcf.gz"), path("${idNormal}.merged.vcf.gz.tbi"), emit: SVsCombinedOutputGermline

  script:
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
    -n -m 1 \\
    -l ${labelparam} \\
    -o ${idNormal}.merged.raw.vcf \\
    -f -d -s -v \\
    ${inVCFs}

  cat ${idNormal}.merged.raw.vcf | \\
    awk -F"\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{print}' | \\
  bcftools sort --temp-dir ./ \\
    > ${idNormal}.merged.clean.anon.vcf

  python ${custom_scripts}/filter-sv-vcf.py \\
    --input ${idNormal}.merged.clean.anon.vcf \\
    --output ${idNormal}.merged.clean.anon.corrected.vcf \\
    --min ${passMin} 

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
    -o ${idNormal}.merged.clean.vcf \\
    ${idNormal}.merged.clean.anon.corrected.vcf

  bcftools view \\
    --samples ${idNormal} \\
    --output-type z \\
    --output-file ${idNormal}.merged.vcf.gz \\
    ${idNormal}.merged.clean.vcf
  
  tabix --preset vcf ${idNormal}.merged.vcf.gz 
  """
}
