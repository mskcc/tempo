process MergeBamsAndMarkDuplicates {
  tag {idSample}

  input:
    tuple val(idSample), path(bam), val(target)

  output:
    tuple val(idSample), path("${idSample}.md.bam"), path("${idSample}.md.bai"), val(target), emit:  mdBams
    path("size.txt"), emit: sizeOutput

  script:

  bamSize = 0
  bam.each{ bamSize = bamSize + it.size()}

  if (workflow.profile == "juno") {
    if(bamSize > 100.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bamSize < 80.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  memMultiplier = params.mem_per_core ? task.cpus : 1
  
  // when increase memory requested from system every time it retries, keep java Xmx steady, in order to give more memory for java garbadge collection
  originalMem = task.attempt ==1 ? task.memory : originalMem
  maxMem = (memMultiplier * originalMem.toString().split(" ")[0].toInteger() - 3)
  maxMem = maxMem < 4 ? 5 : maxMem
  javaOptions = "--java-options '-Xms4000m -Xmx" + maxMem + "g'"
  """
  samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
  gatk MarkDuplicates \
    ${javaOptions} \
    --TMP_DIR ./ \
    --MAX_RECORDS_IN_RAM 50000 \
    --INPUT ${idSample}.merged.bam \
    --METRICS_FILE ${idSample}.bam.metrics \
    --ASSUME_SORT_ORDER coordinate \
    --CREATE_INDEX true \
    --OUTPUT ${idSample}.md.bam

  echo -e "${idSample}\t`du -hs ${idSample}.md.bam`" > size.txt
  """
}
