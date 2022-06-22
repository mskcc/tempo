process SomaticDellyCall {
  tag {idTumor + "__" + idNormal + '@' + svType}

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  
  input:
    tuple val(idTumor), val(idNormal), val(target), val(inputVcfs), val(inputTbis)
    
  output:
    tuple val(idTumor), val(idNormal), val(target), path("${idTumor}__${idNormal}.delly.vcf.gz"), path("${idTumor}__${idNormal}.delly.vcf.gz.tbi")
    
  script:
  """
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${outputPrefix}.delly.vcf.gz \\
    *.delly.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.vcf.gz 
  """
}
