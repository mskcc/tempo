process GermlineMergeSVs {
  tag "${idNormal}"

  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/", mode: params.publishDirMode, pattern: "*.merged.vcf.{gz,gz.tbi}"

  input:
    tuple val(idNormal), val(target), 
      path(dellyVcfs), path(dellyVcfsIndex), 
      path(mantaFile), path(mantaIndex)
    path(custom_scripts)

  output:
    tuple val(idNormal), path("${idNormal}.merged.vcf.gz"), path("${idNormal}.merged.vcf.gz.tbi"), emit: SVsCombinedOutputGermline

  script:
  labelparam = "delly,manta" 
  inVCFs = "${idNormal}.delly.vcf.gz ${mantaFile} "
  """ 
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${idNormal}.delly.vcf.gz \\
    *.delly.vcf.gz 
  tabix --preset vcf ${idNormal}.delly.vcf.gz 

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
    --min 1 

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
