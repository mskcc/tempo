process GermlineMergeDellyAndManta {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/combined_svs/", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    tuple val(idNormal), val(target), path(dellyVcf), path(dellyVcfIndex), path(mantaVcf), path(mantaVcfIndex)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)

  output:
    tuple val("placeHolder"), val("noTumor"), val(idNormal), path("${idNormal}.delly.manta.vcf.gz*"), emit: dellyMantaCombinedOutputGermline
    tuple val("placeHolder"), val("noTumor"), val(idNormal), path("${idNormal}.delly.manta.vcf.gz"), path("${idNormal}.delly.manta.vcf.gz.tbi"), emit: sv4AggregateGermline

  script:
  """ 
  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${idNormal}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz ${idNormal}.manta.vcf.gz

  tabix --preset vcf ${idNormal}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.delly.manta.vcf.gz \
    ${idNormal}.delly.manta.unfiltered.vcf.gz 
    
  tabix --preset vcf ${idNormal}.delly.manta.vcf.gz
  """
}
