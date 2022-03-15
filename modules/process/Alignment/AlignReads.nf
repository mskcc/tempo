process AlignReads {
  tag {fileID + "@" + lane}   // The tag directive allows you to associate each process executions with a custom label

  publishDir "${params.outDir}/bams/${idSample}/fastp", mode: params.publishDirMode, pattern: "*.{html,json}"

  input:
    tuple val(idSample), val(target), path(fastqFile1), val(sizeFastqFile1), path(fastqFile2), val(sizeFastqFile2), val(fileID), val(lane)
    tuple path(genomeFile), path(bwaIndex)

  output:
    tuple val(idSample), path("*.html"), emit: fastPHtml
    tuple val(idSample), path("*.json"), val(fileID), emit: fastPJson4MultiQC
    path("file-size.txt"), emit: laneSize
    tuple val(idSample), val(target), path("*.sorted.bam"), val(fileID), val(lane), path("*.readId"), emit: sortedBam

  script:
  // LSF resource allocation for juno
  // if running on juno, check the total size of the FASTQ pairs in order to allocate the runtime limit for the job, via LSF `bsub -W`
  // if total size of the FASTQ pairs is over 20 GB, use params.maxWallTimeours
  // if total size of the FASTQ pairs is under 12 GB, use 3h. If there is a 140 error, try again with 6h. If 6h doesn't work, try 500h.
  inputSize = sizeFastqFile1 + sizeFastqFile2
  if (workflow.profile == "juno") {
    if (inputSize > 18.GB) {
      task.time = { params.maxWallTime }
    }
    else if (inputSize < 9.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }
  
  // mem --- total size of the FASTQ pairs in MB (max memory `samtools sort` can take advantage of)
  // memDivider --- If mem_per_core is true, use 1. Else, use task.cpus
  // memMultiplier --- If mem_per_core is false, use 1. Else, use task.cpus
  // originalMem -- If this is the first attempt, use task.memory. Else, use `originalMem`
  mem = (inputSize/1024**2).round()
  memDivider = params.mem_per_core ? 1 : task.cpus
  memMultiplier = params.mem_per_core ? task.cpus : 1
  originalMem = task.attempt ==1 ? task.memory : originalMem

  if ( mem < 6 * 1024 / task.cpus ) {
  // minimum total task memory requirment is 6GB because `bwa mem` need this much to run, and increase by 10% everytime retry
      task.memory = { (6 / memMultiplier * (0.9 + 0.1 * task.attempt) + 0.5).round() + " GB" }
      mem = (5.4 * 1024 / task.cpus).round()
  }
  else if ( mem / memDivider * (1 + 0.1 * task.attempt) > originalMem.toMega() ) {
  // if file size is too big, use task.memory as the max mem for this task, and decrease -M for `samtools sort` by 10% everytime retry
      mem = (originalMem.toMega() / memDivider * (1 - 0.1 * task.attempt) + 0.5).round()
  }
  else {
  // normal situation, `samtools sort` -M = inputSize * 2, task.memory is 110% of `samtools sort` and increase by 10% everytime retry
      task.memory = { (mem * memDivider * (1 + 0.1 * task.attempt) / 1024 + 0.5).round() + " GB" }
      mem = mem
  }

  task.memory = task.memory.toGiga() < 1 ? { 1.GB } : task.memory

  filePartNo = fastqFile1.getSimpleName().split("_R1")[-1]
  """
  rgID=`zcat $fastqFile1 | head -1 | tr ':/' '@' | cut -d '@' -f2-5`
  readGroup="@RG\\tID:\${rgID}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
  touch `zcat $fastqFile1 | head -1 | tr ':/\t ' '@' | cut -d '@' -f2-`.readId
  set -e
  set -o pipefail

  fastq1=${fastqFile1}
  fastq2=${fastqFile2}
  if ${params.anonymizeFQ}; then
    ln -s ${fastqFile1} ${idSample}@\${rgID}@R1${filePartNo}.fastq.gz
    ln -s ${fastqFile2} ${idSample}@\${rgID}@R2${filePartNo}.fastq.gz
    fastq1=`echo ${idSample}@\${rgID}@R1${filePartNo}.fastq.gz`
    fastq2=`echo ${idSample}@\${rgID}@R2${filePartNo}.fastq.gz`
  fi

  fastp --html ${idSample}@\${rgID}${filePartNo}.fastp.html --json ${idSample}@\${rgID}${filePartNo}.fastp.json --in1 \${fastq1} --in2 \${fastq2}
  bwa mem -R \"\${readGroup}\" -t ${task.cpus} -M ${genomeFile} \${fastq1} \${fastq2} | samtools view -Sb - > ${idSample}@\${rgID}${filePartNo}.bam

  samtools sort -m ${mem}M -@ ${task.cpus} -o ${idSample}@\${rgID}${filePartNo}.sorted.bam ${idSample}@\${rgID}${filePartNo}.bam
  echo -e "${fileID}@${lane}\t${inputSize}" > file-size.txt
  """
}