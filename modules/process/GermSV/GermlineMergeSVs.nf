process GermlineMergeSVs {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    tuple val(idNormal), val(target), path(dellyVcf), path(dellyVcfIndex), path(mantaVcf), path(mantaVcfIndex)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)

  output:
    tuple val("placeHolder"), val("noTumor"), val(idNormal), path("${idNormal}.delly.manta.vcf.gz"), path("${idNormal}.delly.manta.vcf.gz.tbi"), emit: SVsCombinedOutputGermline

  script:
  """ 
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${idNormal}.delly.vcf.gz \\
    *.delly.vcf.gz 
  tabix --preset vcf ${idNormal}.delly.vcf.gz 

  mergesvvcf \\
    -n -m 1 \\
    -l delly,manta \\
    -o ${idNormal}.delly.manta.raw.vcf \\
    -f -d -s -v \\
    ${idNormal}.delly.vcf.gz ${mantaVcf}

  cat ${idNormal}.delly.manta.raw.vcf | \\
  awk -F"\\t" '\$1 ~ /^#/ && \$1 !~ /^##/ && \$1 !~ /^#CHROM/{next;}{print}' | \\
  bcftools sort --temp-dir ./ \\
    > ${idNormal}.delly.manta.clean.anon.vcf

  python ${custom_scripts}/filter-sv-vcf.py \\
    --input ${outputPrefix}.delly.manta.clean.anon.vcf \\
    --output ${outputPrefix}.delly.manta.clean.anon.corrected.vcf \\
    --min 1 

  bcftools annotate \\
    --set-id 'TEMPO_%INFO/SVTYPE\\_%CHROM\\_%POS' \\
    -o ${outputPrefix}.delly.manta.clean.vcf \\
    ${outputPrefix}.delly.manta.clean.anon.corrected.vcf

  bcftools filter \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.delly.manta.vcf.gz \
    ${idNormal}.delly.manta.unfiltered.vcf.gz 
    
  tabix --preset vcf ${idNormal}.delly.manta.vcf.gz
  """
}
