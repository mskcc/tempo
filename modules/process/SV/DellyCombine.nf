process DellyCombine {
  tag { outputPrefix }

  publishDir "${params.outDir}/${mode}/${outputPrefix}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  
  input:
    tuple val(idTumor), val(idNormal), val(target), path(inputVcfs), path(inputTbis)
    val(mode)
    
  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.delly.vcf.gz"), path("${outputPrefix}.delly.vcf.gz.tbi")
    
  script:
  outputPrefix = [idTumor,idNormal].unique()
  outputPrefix.remove("")
  outputPrefix = outputPrefix.join("__")
  """
  bcftools concat \\
    --allow-overlaps \\
    --output-type z \\
    --output ${outputPrefix}.delly.vcf.gz \\
    *.delly.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.vcf.gz 
  """
}
