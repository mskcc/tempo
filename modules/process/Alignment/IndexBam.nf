process IndexBam {
  tag "${idSample}"

  input:
    tuple val(idSample), path(bam), val(target)

  output:
    tuple val(idSample), path("${bam}"), path("${bam}.bai"), val(target), emit: indexedBam

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  """
}
