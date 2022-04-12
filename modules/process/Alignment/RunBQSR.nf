 process RunBQSR {
    tag {idSample}
    
    publishDir "${params.outDir}/bams/${idSample}", mode: params.publishDirMode, pattern: "*.bam*"

    input:
      tuple val(idSample), path(bam), path(bai), val(target)
      tuple path(genomeFile), path(genomeIndex), path(genomeDict), path(dbsnp), path(dbsnpIndex), path(knownIndels), path(knownIndelsIndex) 

    output:
      tuple val(idSample), val(target), path("${idSample}.bam"), path("${idSample}.bam.bai"), emit: bamsBQSR
      tuple val(idSample), val(target), val("${params.outDir}/bams/${idSample}/${idSample}.bam"), val("${params.outDir}/bams/${idSample}/${idSample}.bam.bai"), emit: bamResults
      path("file-size.txt"), emit: bamSize

    script:
    if (workflow.profile == "juno") {
      if(bam.size() > 100.GB) {
        task.time = { params.maxWallTime }
      }
      else if (bam.size() < 50.GB) {
        task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
      }
      else {
        task.time = task.exitStatus.toString() in params.wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
      }
      task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
    }
    if (task.attempt < 3 ) {
      sparkConf = "Spark --conf 'spark.executor.cores = " + task.cpus + "'"
    }
    else {
      sparkConf=""
      task.cpus = 4
      task.memory = { 6.GB }
      if (workflow.profile == "juno"){ task.time = { params.maxWallTime } }
    }
    
    memMultiplier = params.mem_per_core ? task.cpus : 1
    // when increase memory requested from system every time it retries, keep java Xmx steady, in order to give more memory for java garbadge collection
    originalMem = task.attempt ==1 ? task.memory : originalMem
    maxMem = (memMultiplier * originalMem.toString().split(" ")[0].toInteger() - 3)
    maxMem = maxMem < 4 ? 5 : maxMem
    javaOptions    = "--java-options '-Xmx" + originalMem.toString().split(" ")[0].toInteger() * memMultiplier + "g'"
    knownSites = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
    if ( task.attempt < 3 )
    """
    gatk \
      BQSRPipeline${sparkConf} \
      -R ${genomeFile} \
      -I ${bam} \
      --known-sites ${dbsnp} \
      ${knownSites} \
      --verbosity INFO \
      --create-output-bam-index true \
      -O ${idSample}.bam 
   
    echo -e "${idSample}\t\$(du -b ${idSample}.bam)" > file-size.txt
   
    if [[ -f ${idSample}.bai ]]; then
      mv ${idSample}.bai ${idSample}.bam.bai
    fi
    """
    else 
    """
    gatk \
      BaseRecalibrator${sparkConf} \
      ${javaOptions} \
      --reference ${genomeFile} \
      --known-sites ${dbsnp} \
      ${knownSites} \
      --verbosity INFO \
      --input ${bam} \
      --output ${idSample}.recal.table
    
    gatk \
      ApplyBQSR${sparkConf} \
      ${javaOptions} \
      --reference ${genomeFile} \
      --create-output-bam-index true \
      --bqsr-recal-file ${idSample}.recal.table \
      --input ${bam} \
      --output ${idSample}.bam

    echo -e "${idSample}\t\$(du -b ${idSample}.bam)" > file-size.txt

    if [[ -f ${idSample}.bai ]]; then
      mv ${idSample}.bai ${idSample}.bam.bai
    fi    
    """
  }
