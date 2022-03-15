process generateBasFile {
  tag { idSample }

  input:
  tuple val(idSample), val(target), path(bam), path(bai) 
  path(genomeFile)
  path(genomeIndex)

  output: 
  tuple val(idSample), val(target), path("*.bas")

  script:
  """
  bam_stats \\
    -i ${bam} \\
    -o ${bam}.bas \\
    -r ${genomeIndex} \\
    -@ ${ task.cpus > 1 ? task.cpus - 1 : task.cpus }
  """
}

