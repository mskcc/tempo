process SomaticMergeDellyAndManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(dellyVcfs), path(dellyVcfsIndex), path(mantaFile), path(mantaIndex), file(svabaFile), file(svabaIndex), file(brassFile), file(brassIndex)
    path(custom_scripts) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.delly.manta.vcf.gz"), path("${outputPrefix}.delly.manta.vcf.gz.tbi"), emit: dellyMantaCombined

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  labelparam = "delly,manta" 
  labelparam = svabaFile.name.endsWith("vcf.gz") ?  labelparam + ",svaba" : labelparam 
  labelparam = brassFile.name.endsWith("vcf.gz") ?  labelparam + ",brass" : labelparam 
  """
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${outputPrefix}.delly.vcf.gz \\
    *.delly.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.vcf.gz 

  mergesvvcf \\
    -n -m 1 \\
    -l ${labelparam}\\
    -o ${outputPrefix}.delly.manta.raw.vcf \\
    -f -d -s -v \\
    ${outputPrefix}.delly.vcf.gz ${mantaFile} ${svabaFile} ${brassFile}

  cat ${outputPrefix}.delly.manta.raw.vcf | \\
    awk -F"\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{print}' | \\
  bcftools sort --temp-dir ./ \\
    > ${outputPrefix}.delly.manta.clean.anon.vcf

  python ${custom_scripts}/filter-sv-vcf.py \\
    --input ${outputPrefix}.delly.manta.clean.anon.vcf \\
    --output ${outputPrefix}.delly.manta.clean.anon.corrected.vcf \\
    --min 1 

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
    -o ${outputPrefix}.delly.manta.clean.vcf \\
    ${outputPrefix}.delly.manta.clean.anon.corrected.vcf

  bcftools view \\
    --samples ${idTumor},${idNormal} \\
    --output-type z \\
    --output-file ${outputPrefix}.delly.manta.vcf.gz \\
    ${outputPrefix}.delly.manta.clean.vcf
  
  tabix --preset vcf ${outputPrefix}.delly.manta.vcf.gz 
  """
}
