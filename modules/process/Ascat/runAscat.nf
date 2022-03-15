process runAscat {
  tag {idTumor + "__" + idNormal}
  label 'ascat' 

  input: 
  tuple val(idTumor), val(idNormal), val(target), path(ascatTar), path(tumorBam), path(tumorBai), path(normalBam), path(normalBai) 
  tuple path(genomeFile), path(genomeIndex), path(snpGcCorrections) 

  output:
  tuple val(idTumor), val(idNormal), val(target), file("ascatResults/*.copynumber.caveman.csv"), emit: caveman 
  tuple val(idTumor), val(idNormal), val(target), file("ascatResults/*.samplestatistics.txt"), emit: samplestatistics


  when: params.assayType == "genome"

  script:
  if (params.genome in ["GRCh37","smallGRCh37","GRCh38"]){
    species = "HUMAN"
    assembly = 37
    if (params.genome in ["GRCh38"]) { assembly = 38 }
  } else { // not sure if run will complete with these params. 
    species = params.genome 
    assembly = params.genome
  }
  """
  export TMPDIR=\$(pwd)/tmp 
  mkdir \$TMPDIR 
  for i in ascat_alleleCount_*.tar.gz ;do
    tar -xzf \$i 
  done

  ascat.pl \\
  -o ./ascatResults \\
  -t ${tumorBam} -n ${normalBam} \\
  -sg ${snpGcCorrections} \\
  -r ${genomeFile} \\
  -q 20 -g L \\
  -rs "${species}" -ra "${assembly}" -pr "WGS" \\
  -c ${task.cpus} \\
  -force 
  """
}
