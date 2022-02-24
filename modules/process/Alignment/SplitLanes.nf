process SplitLanesR1 {
  tag {idSample + "@" + fileID}   // The tag directive allows you to associate each process executions with a custom label

  input:
    tuple val(idSample), val(target), file(fastqFile1), val(fileID)

  output:
    path("file-size.txt")
    tuple val(idSample), val(target), path("*R1*.splitLanes.fastq.gz"), path("*.fcid"), path("*.laneCount"), emit: R1SplitData

  when: params.splitLanes

  script:
  inputSize = fastqFile1.size()
  if (workflow.profile == "juno") {
    if (inputSize > 10.GB) {
      task.time = { params.maxWallTime }
    }
    else if (inputSize < 5.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    // if it's the last time to try, use 500h as time limit no matter for what reason it failed before
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  filePartNo = fastqFile1.getSimpleName().split("_R1")[-1]
  filePrefix = fastqFile1.getSimpleName().split("_R1")[0..-2].join("_R1")
  """
    fcid=`zcat $fastqFile1 | head -1 | tr ':/' '@' | cut -d '@' -f2-4`
    touch \${fcid}.fcid
    echo -e "${idSample}@${fileID}\t${inputSize}" > file-size.txt
    zcat $fastqFile1 | awk -v var="\${fcid}" 'BEGIN {FS = ":"} {lane=\$4 ; print | "gzip > ${filePrefix}@"var"_L00"lane"_R1${filePartNo}.splitLanes.fastq.gz" ; for (i = 1; i <= 3; i++) {getline ; print | "gzip > ${filePrefix}@"var"_L00"lane"_R1${filePartNo}.splitLanes.fastq.gz"}}'
    touch `ls *R1*.splitLanes.fastq.gz | wc -l`.laneCount
  """
}

process SplitLanesR2 {
  tag {idSample + "@" + fileID}   // The tag directive allows you to associate each process executions with a custom label

  input:
    tuple val(idSample), val(target), file(fastqFile2), val(fileID)

  output:
    file("file-size.txt")
    tuple val(idSample), val(target), file("*_R2*.splitLanes.fastq.gz"), file("*.fcid"), file("*.laneCount"), emit: R2SplitData

  when: params.splitLanes

  script:
  inputSize = fastqFile2.size()
  if (workflow.profile == "juno") {
    if (inputSize > 10.GB) {
      task.time = { params.maxWallTime }
    }
    else if (inputSize < 5.GB) {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  filePartNo = fastqFile2.getSimpleName().split("_R2")[-1]
  filePrefix = fastqFile2.getSimpleName().split("_R2")[0..-2].join("_R2")
  """
    fcid=`zcat $fastqFile2 | head -1 | tr ':/' '@' | cut -d '@' -f2-4`
    touch \${fcid}.fcid
    echo -e "${idSample}@${fileID}\t${inputSize}" > file-size.txt
    zcat $fastqFile2 | awk -v var="\${fcid}" 'BEGIN {FS = ":"} {lane=\$4 ; print | "gzip > ${filePrefix}@"var"_L00"lane"_R2${filePartNo}.splitLanes.fastq.gz" ; for (i = 1; i <= 3; i++) {getline ; print | "gzip > ${filePrefix}@"var"_L00"lane"_R2${filePartNo}.splitLanes.fastq.gz"}}'
    touch `ls *_R2*.splitLanes.fastq.gz | wc -l`.laneCount
  """
}
