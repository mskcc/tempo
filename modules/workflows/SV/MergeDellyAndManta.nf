process MergeDellyAndManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(dellyVcfs), path(dellyVcfsIndex), path(mantaFile)

  output:
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("${outputPrefix}.delly.manta.vcf.gz"), path("${outputPrefix}.delly.manta.vcf.gz.tbi"), emit: sv4Aggregate

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """ 
  bcftools view \
    --samples ${idTumor},${idNormal} \
    --output-type z \
    --output-file ${outputPrefix}.manta.swap.vcf.gz \
    ${outputPrefix}.manta.vcf.gz 
    
  tabix --preset vcf ${outputPrefix}.manta.swap.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${outputPrefix}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz ${outputPrefix}.manta.swap.vcf.gz
  
  tabix --preset vcf ${outputPrefix}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    ${outputPrefix}.delly.manta.unfiltered.vcf.gz | \
  bcftools sort \
    --output-type z \
    --output-file ${outputPrefix}.delly.manta.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.manta.vcf.gz 
  """
}
