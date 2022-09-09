process CreateBaitsetFiles {
  tag "${targetId}"

  input:
    tuple val(targetId), path("raw_targets.bed"), path("raw_baits.bed")
    path(genomeFile)
    path(genomeIndex)
    path(genomeDict)
    path(codingRegions)
    
  output:
    tuple val(targetId), path(targetInterval), path(baitInterval), emit: baitsetInterval
    tuple val(targetId), path(codingBaitsetBed), emit: codingBaitsetBed
    tuple val(targetId), path("${targetPlus5}.gz"), path("${targetPlus5}.gz.tbi"), emit:baitsetPlus5
    tuple val(targetId), path(targetPlus5), emit:baitsetPlus5_unzipped

  script:
  targetInterval = "${targetId}.targets.ilist"
  baitInterval = "${targetId}.baits.ilist"
  codingBaitsetBed = "${targetId}.coding.bed"
  targetBed = "${targetId}.targets.bed"
  baitBed = "${targetId}.baits.bed"
  targetPlus5 = "${targetId}.plus5bp.bed"
  """
  bedtools sort -i raw_targets.bed | bedtools merge -i - > ${targetBed}
  bedtools sort -i raw_baits.bed | bedtools merge -i - > ${baitBed}

  bedtools intersect \\
    -a ${codingRegions} \\
    -b ${targetBed} > \\
    intersect.bed 
  sort -k1,1 -k 2,2n -k 3,3n intersect.bed > intersect.sorted.bed
  bedtools merge -i intersect.sorted.bed > ${codingBaitsetBed} 

  cut -f 1,2 ${genomeIndex} > this.genome
  bedtools slop \\
    -i ${targetBed} \\
    -g ./this.genome \\
    -b 5 > \\
    ${targetPlus5}
  bgzip -c ${targetPlus5} > ${targetPlus5}.gz
  tabix -p bed ${targetPlus5}.gz

  gatk BedToIntervalList \\
    -I ${targetBed} \\
    -O ${targetInterval} \\
    -SD ${genomeDict}

  gatk BedToIntervalList \\
    -I ${baitBed} \\
    -O ${baitInterval} \\
    -SD ${genomeDict}

  """
}
