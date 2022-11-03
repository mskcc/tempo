process GermlineMergeSVs {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.merged.vcf.{gz,gz.tbi}"
  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/intermediate_files", mode: params.publishDirMode, pattern: "*.merged.raw.vcf.{gz,gz.tbi}"

  input:
    tuple val(idNormal), val(target), 
      path(Vcfs), path(Tbis),
      val(callerNames)
    path(custom_scripts)

  output:
    tuple val(idNormal), val(target), path("${idNormal}.merged.vcf.gz"), path("${idNormal}.merged.vcf.gz.tbi"), emit: SVsCombinedOutputGermline
    path("${idNormal}.merged.raw.vcf.{gz,gz.tbi}")

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
    awk -F"\\t" -v OFS="\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = "."; print \$0}' | \\
  bcftools sort --temp-dir ./ \\
    > ${idNormal}.merged.clean.anon.vcf

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
    -o ${idNormal}.merged.clean.vcf \\
    ${idNormal}.merged.clean.anon.vcf

  python ${custom_scripts}/filter-sv-vcf.py \\
    --input ${idNormal}.merged.clean.vcf \\
    --output ${idNormal}.merged.clean.corrected.vcf \\
    --min ${passMin}

  bcftools view \\
    --samples ${idNormal} \\
    --output-type z \\
    --output-file ${idNormal}.merged.vcf.gz \\
    ${idNormal}.merged.clean.corrected.vcf
  
  tabix --preset vcf ${idNormal}.merged.vcf.gz

  bcftools view -O z -o ${idNormal}.merged.raw.vcf.gz ${idNormal}.merged.raw.vcf
  """
}
