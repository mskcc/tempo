process runAscatAlleleCount {
  tag {idTumor + "__" + idNormal + "@" + ascatIndex }
  label 'ascat' 

  input: 
  each ascatIndex 
  val(ascatIndexLimit)
  tuple val(idTumor), val(idNormal), val(target), path(tumorBam), path(tumorBai), path(normalBam), path(normalBai) 
  path(genomeFile)
  path(genomeIndex)
  path(snpGcCorrections) 
  
  output:
  tuple val(idTumor), val(idNormal), val(target), path("ascat_alleleCount_${ascatIndex}.tar.gz")

  when: params.assayType == "genome"

  script:
  genome = params.genome
  if (genome in ["GRCh37","smallGRCh37","GRCh38"]){
    species = "HUMAN"
    assembly = 37
    if (genome in ["GRCh38"]) { assembly = 38 }
  } else { // not sure if run will complete with these params. 
    species = genome 
    assembly = genome
  }
  if (params.ascatAlleleCount == 1 ) { 
    indexParam = ""
  } else {
    indexParam = "-i ${ascatIndex} -x ${params.ascatAlleleCount}"
  }
  """
  mkdir -p ascatResults
  export TMPDIR=\$(pwd)/tmp 
  mkdir \$TMPDIR 

  ascat.pl \\
  -o ./ascatResults \\
  -t ${tumorBam} -n ${normalBam} \\
  -sg ${snpGcCorrections} \\
  -r ${genomeFile} \\
  -q 20 -g L \\
  -rs "${species}" -ra "${assembly}" -pr "WGS" \\
  -c ${task.cpus} \\
  -force \\
  -p allele_count ${indexParam}
  
  tar -czf ascat_alleleCount_${ascatIndex}.tar.gz ascatResults/
  """
}
