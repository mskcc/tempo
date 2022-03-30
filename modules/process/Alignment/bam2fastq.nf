process bam2fastq {
  tag "${idSample}"

  input:
    tuple val(idSample), val(target), file(bam), val(bai)

  output:
    tuple val(idSample), val(target), path("*_1.fastq.gz"), path("*_2.fastq.gz"), emit: fastqOutput

  script:
  inputSize = bam.size()
  if (workflow.profile == "juno") {
    if (inputSize > 80.GB) {
      task.time = { params.maxWallTime }
    }
    else if (inputSize < 40.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    // if it's the last time to try, use 500h as time limit no matter for what reason it failed before
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }
    memMultiplier = params.mem_per_core ? task.cpus : 1
    // when increase memory requested from system every time it retries, keep java Xmx steady, in order to give more memory for java garbadge collection
    originalMem = task.attempt ==1 ? task.memory : originalMem
    maxMem = (memMultiplier * originalMem.toString().split(" ")[0].toInteger() - 3)
    maxMem = maxMem < 4 ? 5 : maxMem
    javaOptions    = "--java-options '-Xmx" + originalMem.toString().split(" ")[0].toInteger() * memMultiplier + "g'"

  """
  gatk SamToFastq ${javaOptions} VALIDATION_STRINGENCY=LENIENT I=${bam} RG_TAG=ID OUTPUT_PER_RG=true COMPRESS_OUTPUTS_PER_RG=true OUTPUT_DIR=./ INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=true
  ls *.fastq.gz | xargs -I {} -n1 mv {} `basename ${idSample}`@{}
  """
}
