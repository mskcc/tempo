#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
================================================================================
--------------------------------------------------------------------------------
Processes overview:

Alignment
----------------
 - AlignReads
    --- Map paired-end FASTQs with bwa mem
    --- Sort BAM with samtools sort
    --- FASTQ QC with FastP 
 - MergeBam --- Merge BAM for the same samples from different fileID, samtools merge
 - MarkDuplicates --- Mark Duplicates with GATK4 MarkDuplicates
 - CreateRecalibrationTable --- Create Recalibration Table with GATK4 BaseRecalibrator
 - RecalibrateBam --- Recalibrate Bam with GATK4 ApplyBQSR

Somatic Analysis
----------------
 - CreateScatteredIntervals --- GATK4 SplitIntervals
 - RunMutect2 --- somatic SNV calling, MuTect2
 - SomaticCombineMutect2Vcf --- combine Mutect2 calls, bcftools
 - SomaticDellyCall --- somatic SV calling, Delly
 - SomaticRunManta --- somatic SV calling, Manta
 - SomaticMergeDellyAndManta --- combine Manta and Delly VCFs
 - SomaticRunStrelka2 --- somatic SNV calling, Strelka2, using Manta for small indel calling by default
 - SomaticCombineChannel --- combine and filter VCFs, bcftools
 - SomaticAnnotateMaf --- annotate MAF, vcf2maf
 - RunMutationSignatures --- mutational signatures
 - DoFacets --- facets-suite: mafAnno.R, geneLevel.R, armLevel.R
 - RunPolysolver --- Polysolver
 - RunLOHHLA --- LOH in HLA
 - SomaticFacetsAnnotation --- annotate FACETS
 - RunNeoantigen --- NetMHCpan 4.0
 - RunMsiSensor --- MSIsensor
 - MetaDataParser --- python script to parse metadata into single *tsv

Germline Analysis
-----------------
 - CreateScatteredIntervals --- (run once) GATK4 SplitIntervals
 - GermlineDellyCall --- germline SV calling and filtering, Delly
 - GermlineRunManta --- germline SV calling, Manta
 - GermlineMergeDellyAndManta --- merge SV calls from Delly and Manta
 - GermlineRunHaplotypecaller --- germline SNV calling, GATK4
 - GermlineCombineHaplotypecallerVcf --- concatenate VCFs of GATK4 HaplotypeCaller
 - GermlineRunStrelka2 --- germline SNV calling, Strelka2 (with InDels from Manta)
 - GermlineCombineChannel --- combined and filter germline calls, bcftools
 - GermlineAnnotateMaf--- annotate MAF, vcf2maf

Quality Control
-----------------
 - QcAlfred - BAM QC metrics
 - QcCollectHsMetrics --- *For WES only* Calculate hybrid-selection metrics, GATK4 CollectHsMetrics
 - QcBamAggregate --- aggregates information from QcAlfred and QcCollectHsMetrics across all samples
 - QcConpair --- Tumor-Normal quality/contamination
 - QcConpairAll --- Tumor-Normal All Combination quality/contamination

Cohort Aggregation
-----------------
 - SomaticAggregateMaf --- collect outputs, MAF
 - SomaticAggregateNetMHC --- collect outputs, neoantigen prediction
 - SomaticAggregateFacets --- collect outputs, FACETS
 - SomaticAggregateSv --- collect outputs, SVs
 - SomaticAggregateLOHHLA --- collect outputs, LOHHLA
 - SomaticAggregateMetaData --- collect outputs, sample data
 - GermlineAggregateMaf --- collect outputs, MAF
 - GermlineAggregateSv --- collect outputs, SVs
 - QcConpairAggregate --- aggregates information from QcConpair or QcConpairAll across all sample

*/

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println "ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)"
  exit 1
}

// User-set runtime parameters
publishAll = params.publishAll
outDir     = file(params.outDir).toAbsolutePath()
outname = params.outname
runGermline = params.germline
runSomatic = params.somatic
runQC = params.QC
runAggregate = params.aggregate
runConpairAll = false
wallTimeExitCode = params.wallTimeExitCode ? params.wallTimeExitCode.split(',').collect{it.trim().toLowerCase()} : []
multiqcWesConfig = workflow.projectDir + "/lib/multiqc_config/exome_multiqc_config.yaml"
multiqcWgsConfig = workflow.projectDir + "/lib/multiqc_config/wgs_multiqc_config.yaml"
multiqcTempoLogo = workflow.projectDir + "/docs/tempoLogo.png"
epochMap = [:]
if (params.watch == true){
  touchInputs()
}

println ""

pairingQc = params.pairing

if (params.mapping || params.bamMapping) {
  TempoUtils.checkAssayType(params.assayType)
  if (params.watch == false) {
    mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)
    (checkMapping1, checkMapping2, inputMapping) = params.mapping ? TempoUtils.extractFastq(mappingFile, params.assayType).into(3) : TempoUtils.extractBAM(mappingFile, params.assayType).into(3)
  }
  else if (params.watch == true) {
    mappingFile = params.mapping ? file(params.mapping, checkIfExists: false) : file(params.bamMapping, checkIfExists: false)
    (checkMapping1, checkMapping2, inputMapping) = params.mapping ? watchMapping(mappingFile, params.assayType).into(3) : watchBamMapping(mappingFile, params.assayType).into(3)
    epochMap[params.mapping ? params.mapping : params.bamMapping ] = 0
  }
  else{}
  if(params.pairing){
    if(runQC){
      pairingQc = true
    }
    if (params.watch == false) {
      pairingFile = file(params.pairing, checkIfExists: true)
      (checkPairing1, checkPairing2, inputPairing) = TempoUtils.extractPairing(pairingFile).into(3)
      TempoUtils.crossValidateTargets(checkMapping1, checkPairing1)
      if(!TempoUtils.crossValidateSamples(checkMapping2, checkPairing2)){exit 1}
    }
    else if (params.watch == true) {
      pairingFile = file(params.pairing, checkIfExists: false)
      (checkPairing1, checkPairing2, inputPairing) = watchPairing(pairingFile).into(3)
      epochMap[params.pairing] = 0 
    }
    else{}
    if (!runSomatic && !runGermline && !runQC){
      println "ERROR: --pairing [tsv] is not used because none of --somatic/--germline/--QC was enabled. If you only need to do BAM QC and/or BAM generation, remove --pairing [tsv]."
      exit 1
    }
  }
  else{
    if (runSomatic || runGermline){
      println "ERROR: --pairing [tsv] needed when using --mapping/--bamMapping [tsv] with --somatic/--germline"
      exit 1
    }
  }
}
else{
  if(params.pairing){
    println "ERROR: When --pairing [tsv], --mapping/--bamMapping [tsv] must be provided."
    exit 1
  }
}

if (!runSomatic && runGermline){
    println "WARNING: You can't run GERMLINE section without running SOMATIC section. Activating SOMATIC section automatically"
    runSomatic = true
}

if (runAggregate == false){
  if (!params.mapping && !params.bamMapping){
    println "ERROR: (--mapping/-bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv] ) or (--aggregate [tsv]) need to be provided, otherwise nothing to be run."
    exit 1
  }
}
else if (runAggregate == true){
  if ((params.mapping || params.bamMapping) && params.pairing){
    if (!(runSomatic || runGermline || runQC)){
      println "ERROR: Nothing to be aggregated. One or more of the option --somatic/--germline/--QC need to be enabled when using --aggregate"
    }
  }
  else if ((params.mapping || params.bamMapping) && !params.pairing){
    if (!runQC){
      println "ERROR: Nothing to be aggregated. --QC need to be enabled when using --mapping/--bamMapping [tsv], --pairing false and --aggregate true."
      exit 1
    }
  }
  else{
    println "ERROR: (--mapping/--bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv]) or (--aggregate [tsv]) need to be provided when using --aggregate true"
    println "       If you want to run aggregate only, you need to use --aggregate [tsv]. See manual"
    exit 1
  }

}
else {
  if ((runSomatic || runGermline || runQC) && !params.mapping && !params.bamMapping){
    println "ERROR: Conflict input! When running --aggregate [tsv] with --mapping/--bamMapping/--pairing [tsv] disabled, --QC/--somatic/--germline all need to be disabled!"
    println "       If you want to run aggregate somatic/germline/qc, just include an additianl colum PATH in the [tsv] and no need to use --QC/--somatic/--germline flag, since it's auto detected. See manual"
    exit 1
  }
}

if (!(params.cosmic in ['v2', 'v3'])) {
  println "ERROR: Possible values of mutational signature reference --cosmic is 'v2', 'v3'"
  exit 1
}

referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


// Skip these processes if starting from aligned BAM files
if (params.mapping) {

  // Parse input FASTQ mapping
  if(params.watch != true){
  inputMapping.groupTuple(by: [0])
              .map { idSample, targets, files_pe1, files_pe2
                -> tuple(groupKey(idSample, targets.size()), targets, files_pe1, files_pe2)
              }
              .transpose()
	      .set{ inputMapping }
  }

  inputMapping.map{ idSample, target, file_pe1, file_pe2 ->
                   [idSample, target, file_pe1, file_pe2, idSample + "@" + file_pe1.getSimpleName(), file_pe1.getSimpleName()]
              }
              .set{ inputFastqs }

  if (params.splitLanes) {
  inputFastqs
        .into{ fastqsNeedSplit; fastqsNoNeedSplit }

  fastqsNeedSplit
        .filter{ item -> !(item[2].getName() =~ /_L(\d){3}_/) }
        .multiMap{ idSample, target, file_pe1, file_pe2, fileID, lane ->
		inputFastqR1: [idSample, target, file_pe1, file_pe1.toString()]
		inputFastqR2: [idSample, target, file_pe2, file_pe2.toString()]
	}
        .set{ fastqsNeedSplit }

  fastqsNoNeedSplit
        .filter{ item -> item[2].getName() =~ /_L(\d){3}_/ }
        .map { idSample, target, file_pe1, file_pe2, fileID, lane
                -> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
        }
        .set{ fastqsNoNeedSplit }


  process SplitLanesR1 {
    tag {idSample + "@" + fileID}   // The tag directive allows you to associate each process executions with a custom label

    input:
      set idSample, target, file(fastqFile1), fileID from fastqsNeedSplit.inputFastqR1

    output:
      file("file-size.txt") into R1Size
      set idSample, target, file("*R1*.splitLanes.fastq.gz"), file("*.fcid"), file("*.laneCount") into perLaneFastqsR1

    when: params.splitLanes

    script:
    inputSize = fastqFile1.size()
    if (workflow.profile == "juno") {
      if (inputSize > 10.GB) {
        task.time = { params.maxWallTime }
      }
      else if (inputSize < 5.GB) {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
      }
      else {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
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
      set idSample, target, file(fastqFile2), fileID from fastqsNeedSplit.inputFastqR2

    output:
      file("file-size.txt") into R2Size
      set idSample, target, file("*_R2*.splitLanes.fastq.gz"), file("*.fcid"), file("*.laneCount") into perLaneFastqsR2

    when: params.splitLanes

    script:
    inputSize = fastqFile2.size()
    if (workflow.profile == "juno") {
      if (inputSize > 10.GB) {
        task.time = { params.maxWallTime }
      }
      else if (inputSize < 5.GB) {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
      }
      else {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
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

  def fastqR1fileIDs = [:]
  perLaneFastqsR1 = perLaneFastqsR1.transpose()
        .map{ item ->
          def idSample = item[0]
	  def target = item[1]
	  def fastq = item[2]
	  def fileID = idSample + "@" + item[3].getSimpleName()
	  def lane = fastq.getSimpleName().split("_L00")[1].split("_")[0]
	  def laneCount = item[4].getSimpleName().toInteger()

	  // This only checks if same read groups appears in two or more fastq files which belongs to the same sample. Cross sample check will be performed after AlignReads since the read group info is not available for fastqs which does not need to be split.
	  if ( !params.watch ){
	  if(!TempoUtils.checkDuplicates(fastqR1fileIDs, fileID + "@" + lane, fileID + "\t" + fastq, "the follwoing fastq files since they contain the same RGID")){exit 1}
	  }

	  [idSample, target, fastq, fileID, lane, laneCount]
        }

  def fastqR2fileIDs = [:]
  perLaneFastqsR2 = perLaneFastqsR2.transpose()
        .map{ item ->
          def idSample = item[0]
	  def target = item[1]
	  def fastq = item[2]
	  def fileID = idSample + "@" + item[3].getSimpleName()
	  def lane = fastq.getSimpleName().split("_L00")[1].split("_")[0]
	  def laneCount = item[4].getSimpleName().toInteger()

	  if ( !params.watch ){
	  if(!TempoUtils.checkDuplicates(fastqR2fileIDs, fileID + "@" + lane, fileID + "\t" + fastq, "the follwoing fastq files since they contain the same RGID")){exit 1}
	  }

	  [idSample, target, fastq, fileID, lane, laneCount]
        }

  fastqFiles  = perLaneFastqsR1
        .mix(perLaneFastqsR2)
        .groupTuple(by: [0,1,3,4,5], size: 2, sort: true)
        .map {  idSample, target, fastqPairs, fileID, lanes, laneCount ->
          tuple(idSample, target, fastqPairs, groupKey(fileID, laneCount), lanes)
        }
        .map{ idSample, target, fastqPairs, fileID, lane ->
             [idSample, target, fastqPairs[0], fastqPairs[1], fileID, lane]
        }
        .map{ item ->
                def idSample = item[0]
                def target = item[1]
                def fastqPair1 = item[2]
                def fastqPair2 = item[3]
                if (item[2].toString().split("_R1").size() < item[3].toString().split("_R1").size()) {
                  fastqPair1 = item[3]
                  fastqPair2 = item[2]
                }
                def fileID = item[4]
                def lane = item[5]

             [idSample, target, fastqPair1, fastqPair1.size(), fastqPair2, fastqPair2.size(), fileID, lane]
        }
        .mix(fastqsNoNeedSplit)
  }
  else{
     fastqFiles =  inputFastqs.map { idSample, target, file_pe1, file_pe2, fileID, lane
					-> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
				   }
  }

  // AlignReads - Map reads with BWA mem output SAM
  process AlignReads {
    tag {fileID + "@" + lane}   // The tag directive allows you to associate each process executions with a custom label

    publishDir "${outDir}/bams/${idSample}/fastp", mode: params.publishDirMode, pattern: "*.{html,json}"

    input:
      set idSample, target, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, fileID, lane from fastqFiles
      set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

    output:
      set idSample, file("*.html") into fastPHtml
      set idSample, file("*.json"), fileID into fastPJson4MultiQC
      file("file-size.txt") into laneSize
      set idSample, target, file("*.sorted.bam"), fileID, lane, file("*.readId") into sortedBam

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
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
      }
      else {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
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

  fastPJson4MultiQC
    .groupTuple(by:[2])
    .map{idSample, jsonFile, fileID -> 
      def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
      [idSampleout, jsonFile]
    }.groupTuple(by: [0])
    .map{ idSample, jsonFile -> 
      [idSample, jsonFile.flatten()]
    }.into{fastPJson4cohortMultiQC; fastPJson4sampleMultiQC} 

  // Check for FASTQ files which might have different path but contains the same reads, based only on the name of the first read.
  def allReadIds = [:]
  sortedBam.map { idSample, target, bam, fileID, lane, readIdFile -> def readId = "@" + readIdFile.getSimpleName().replaceAll("@", ":")

		// Use the first line of the fastq file (the name of the first read) as unique identifier to check across all the samples if there is any two fastq files contains the same read name, if so, we consider there are some human error of mixing up the same reads into different fastq files
		if ( !params.watch ){
		if(!TempoUtils.checkDuplicates(allReadIds, readId, idSample + "\t" + bam, "the follwoing samples, since they contain the same read: \n${readId}")){exit 1}
		}

		[idSample, target, bam, fileID, lane]
	   }
	   .groupTuple(by: [3])
	   .map{ item ->
		      def idSample = item[0] instanceof Collection ? item[0].first() : item[0]
		      def target   = item[1] instanceof Collection ? item[1].first() : item[1]
		      def bams = item[2]
		      [idSample, target, bams]
	   }
	   .groupTuple(by: [0])
	   .map{ item ->
		      def idSample = item[0]
		      def target =  item[1] instanceof Collection ? item[1].first() : item[1]
		      def bams = item[2].flatten()
		      [idSample, bams, target]
	   }
	   .set{ groupedBam }

  // MergeBams
  process MergeBams {
    tag {idSample}

    input:
      set idSample, file(bam), target from groupedBam

    output:
      set idSample, file("${idSample}.merged.bam"), target into mergedBam
      file("size.txt") into sizeOutput

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
    echo -e "${idSample}\t`du -hs ${idSample}.merged.bam`" > size.txt
    """
  }


// GATK MarkDuplicates
  process MarkDuplicates {
    tag {idSample}

    input:
      set idSample, file(bam), target from mergedBam

    output:
      set idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), target into mdBams, mdBams4BQSR

    script:
    if (workflow.profile == "juno") {
      if(bam.size() > 120.GB) {
        task.time = { params.maxWallTime }
      }
      else if (bam.size() < 100.GB) {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
      }
      else {
        task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
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
    gatk MarkDuplicates \
      ${javaOptions} \
      --TMP_DIR ./ \
      --MAX_RECORDS_IN_RAM 50000 \
      --INPUT ${idSample}.merged.bam \
      --METRICS_FILE ${idSample}.bam.metrics \
      --ASSUME_SORT_ORDER coordinate \
      --CREATE_INDEX true \
      --OUTPUT ${idSample}.md.bam
    """
  }


 // GATK BaseRecalibrator , CreateRecalibrationTable 
  process CreateRecalibrationTable {
    tag {idSample}

    input:
      set idSample, file(bam), file(bai), target from mdBams
      set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex) from Channel.value([
        referenceMap.genomeFile,
        referenceMap.genomeIndex,
        referenceMap.genomeDict,
        referenceMap.dbsnp,
        referenceMap.dbsnpIndex,
        referenceMap.knownIndels,
        referenceMap.knownIndelsIndex 
      ])

    output:
      set idSample, file("${idSample}.recal.table") into recalibrationTable

    script:
    if (task.attempt < 3 ){
      sparkConf = " BaseRecalibratorSpark --conf 'spark.executor.cores = " + task.cpus + "'"
      if (workflow.profile == "juno") {
        if (bam.size() > 480.GB) {
          task.time = { params.maxWallTime }
        }
        else if (bam.size() < 240.GB) {
          task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
        }
        else {
          task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
        }
      }
    }
    else {
      sparkConf = " BaseRecalibrator"
      task.cpus = 4
      task.memory = { 6.GB }
      if (workflow.profile == "juno"){ task.time = { params.maxWallTime } }
    }

    memMultiplier = params.mem_per_core ? task.cpus : 1
    // when increase memory requested from system every time it retries, keep java Xmx steady, in order to give more memory for java garbadge collection
    originalMem = task.attempt ==1 ? task.memory : originalMem
    javaOptions = "--java-options '-Xmx" + originalMem.toString().split(" ")[0].toInteger() * memMultiplier + "g'"

    knownSites = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
    """
    gatk \
      ${sparkConf} \
      ${javaOptions} \
      --reference ${genomeFile} \
      --known-sites ${dbsnp} \
      ${knownSites} \
      --verbosity INFO \
      --input ${bam} \
      --output ${idSample}.recal.table
    """ 
  }

  mdBams4BQSR.combine(recalibrationTable, by:[0]).set{ inputsBQSR }

  // GATK ApplyBQSR, RecalibrateBAM
  process RecalibrateBam {
    tag {idSample}

    publishDir "${outDir}/bams/${idSample}", mode: params.publishDirMode, pattern: "*.bam*"

    input:
      set idSample, file(bam), file(bai), target, file(recalibrationReport) from inputsBQSR
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict 
      ])

    output:
      set idSample, target, file("${idSample}.bam"), file("${idSample}.bam.bai") into bamsBQSR4Alfred, bamsBQSR4CollectHsMetrics, bamsBQSR4Tumor, bamsBQSR4Normal, bamsBQSR4QcPileup, bamsBQSR4Qualimap
      set idSample, target, val("${file(outDir).toString()}/bams/${idSample}/${idSample}.bam"), val("${file(outDir).toString()}/bams/${idSample}/${idSample}.bam.bai") into bamResults
      file("file-size.txt") into bamSize

    script:

    if (task.attempt < 3 ) {
      sparkConf = " ApplyBQSRSpark --conf 'spark.executor.cores = " + task.cpus + "'"
      if (workflow.profile == "juno") {
        if (bam.size() > 200.GB){
          task.time = { params.maxWallTime }
        }
        else if (bam.size() < 100.GB) {
          task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
        }
        else {
          task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
        }
      }
    }
    else {
      sparkConf = " ApplyBQSR"
      task.cpus = 4
      task.memory = { 6.GB }
      if (workflow.profile == "juno"){ task.time = { params.maxWallTime } }
    }

    memMultiplier = params.mem_per_core ? task.cpus : 1
    // when increase memory requested from system every time it retries, keep java Xmx steady, in order to give more memory for java garbadge collection
    originalMem = task.attempt ==1 ? task.memory : originalMem
    javaOptions = "--java-options '-Xmx" + originalMem.toString().split(" ")[0].toInteger() * memMultiplier + "g'"
    """
    echo -e "${idSample}\t${bam.size()}" > file-size.txt
    gatk \
      ${sparkConf} \
      ${javaOptions} \
      --reference ${genomeFile} \
      --create-output-bam-index true \
      --bqsr-recal-file ${recalibrationReport} \
      --input ${bam} \
      --output ${idSample}.bam
    if [[ -f ${idSample}.bai ]]; then
      mv ${idSample}.bai ${idSample}.bam.bai
    fi
    """
  }

  File file = new File(outname)
  file.newWriter().withWriter { w ->
      w << "SAMPLE\tTARGET\tBAM\tBAI\n"
  }

  bamResults.subscribe { Object obj ->
      file.withWriterAppend { out ->
          out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}"
      }
  }
} // End of "if (params.mapping) {}"


/*
================================================================================
=                                PAIRING TUMOR and NORMAL                      =
================================================================================
*/

// If starting with BAM files, parse BAM pairing input
if (params.bamMapping) {
  inputMapping.into{bamsBQSR4Alfred; bamsBQSR4Qualimap; bamsBQSR4CollectHsMetrics; bamsBQSR4Tumor; bamsBQSR4Normal; bamsBQSR4QcPileup; bamPaths4MultiQC}

  if (runQC){
    bamPaths4MultiQC.map{idSample, target, bam, bai ->
      [ idSample,target, bam.getParent() ]
    }.set{ locateFastP4MultiQC}

    locateFastP4MultiQC.map{ idSample,target, bamFolder -> 
      [idSample, file(bamFolder + "/fastp/*json")]
    }.into{fastPJson4sampleMultiQC;fastPJson4cohortMultiQC}
  } 
}

if (params.pairing) {

  // Parse input FASTQ mapping

  inputPairing.into{pairing4T; pairing4N; pairingTN; pairing4QC; inputPairing}
  bamsBQSR4Tumor.combine(pairing4T)
                          .filter { item ->
                            def idSample = item[0]
                            def target = item[1]
                            def sampleBam = item[2]
                            def sampleBai = item[3]
                            def idTumor = item[4]
                            def idNormal = item[5]
                            idSample == idTumor
                          }.map { item ->
                            def idTumor = item[4]
                            def idNormal = item[5]
                            def tumorBam = item[2]
                            def tumorBai = item[3]
                            def target = item[1]
                            return [ idTumor, idNormal, target, tumorBam, tumorBai ]
                          }
			  .unique()
			  .into{bamsTumor4Combine; bamsTumor4VcfCombine}

  bamsBQSR4Normal.combine(pairing4N)
                          .filter { item ->
                            def idSample = item[0]
                            def target = item[1]
                            def sampleBam = item[2]
                            def sampleBai = item[3]
                            def idTumor = item[4]
                            def idNormal = item[5]
                            idSample == idNormal
                          }.map { item ->
                            def idTumor = item[4]
                            def idNormal = item[5]
                            def normalBam = item[2]
                            def normalBai = item[3]
                            def target = item[1]
                            return [ idTumor, idNormal, target, normalBam, normalBai ]
                          }.unique()
			  .into{ bamsNormal4Combine; bamsNormalOnly }
  bamsNormalOnly.map { item ->
			def idNormal = item[1]
			def target = item[2]
			def normalBam = item[3]
			def normalBai = item[4]
			return [ idNormal, target, normalBam, normalBai ] }
	.unique()
	.into{ bams4Haplotypecaller; bamsNormal4Polysolver; bamsForStrelkaGermline; bamsForMantaGermline; bamsForDellyGermline }


  bamsTumor4Combine.combine(bamsNormal4Combine, by: [0,1,2])
                          .map { item -> // re-order the elements
                            def idTumor = item[0]
                            def idNormal = item[1]
                            def target = item[2]
                            def bamTumor = item[3]
                            def baiTumor = item[4]
                            def bamNormal = item[5]
                            def baiNormal = item[6]

                            return [ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal ]
                          }
			  .set{bamFiles}


} // End of "if (pairingQc) {}"

/*
================================================================================
=                                SOMATIC PIPELINE                              =
================================================================================
*/

if (runSomatic || runGermline || runQC) {
  // parse --tools parameter for downstream 'when' conditionals, e.g. when: 'delly ' in tools
  tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

  // Allow shorter names
  if ("mutect" in tools) {
    tools.add("mutect2")
  }
  if ("strelka" in tools) {
    tools.add("strelka2")
  }

  // If using Strelka2, run Manta as well to generate candidate indels
  if ("strelka2" in tools) {
    tools.add("manta")
  }

  // If using running either conpair or conpairAll, run pileup as well to generate pileups
  if ("conpair" in tools) {
    tools.add("pileup")
  }
  if ("conpairall" in tools) {
    runConpairAll = true
    tools.add("pileup")
  }

}

if (runSomatic || runGermline) {
// GATK SplitIntervals, CreateScatteredIntervals
process CreateScatteredIntervals {

  input:
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
      ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
      ])
  
  output:
    set file("agilent*.interval_list"), val("agilent"), val("agilent") into agilentIList
    set file("idt*.interval_list"), val("idt"), val("idt") into idtIList
    set file("wgs*.interval_list"), val("wgs"), val("wgs") into wgsIList

  when: runSomatic || runGermline

  script:
  scatterCount = params.scatterCount
  """
  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${agilentTargets} \
    --scatter-count ${scatterCount} \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output agilent

  for i in agilent/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i agilent-\$BASENAME
  done

  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${idtTargets} \
    --scatter-count ${scatterCount} \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output idt

  for i in idt/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i idt-\$BASENAME
  done

  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${wgsTargets} \
    --scatter-count ${scatterCount} \
    --subdivision-mode INTERVAL_SUBDIVISION \
    --output wgs 

  for i in wgs/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i wgs-\$BASENAME
  done
  """
}


agilentIList.mix(idtIList, wgsIList).into{mergedIList4T; mergedIList4N}

//Associating interval_list files with BAM files, putting them into one channel

bamFiles.into{bamsTN4Intervals; bamsForDelly; bamsForManta; bams4Strelka; bamns4CombineChannel; bamsForMsiSensor; bamFiles4DoFacets; bamsForLOHHLA }

bamsTN4Intervals.combine(mergedIList4T, by: 2).map{
  item ->
    def idTumor = item[1]
    def idNormal = item[2]
    def target = item[0]
    def tumorBam = item[3]
    def normalBam = item[4]
    def tumorBai = item[5]
    def normalBai = item[6]
    def intervalBed = item[7]
    def key = idTumor+"__"+idNormal+"@"+target // adding one unique key

    return [ key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed ]
}.map{ 
    key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed -> 
    tuple ( 
         groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
         idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed
    )
}
.transpose()
.set{ mergedChannelSomatic }


bams4Haplotypecaller.combine(mergedIList4N, by: 1)
.map{
  item ->
    def idNormal = item[1]
    def target = item[0]
    def normalBam = item[2]
    def normalBai = item[3]
    def intervalBed = item[4]
    def key = idNormal+"@"+target // adding one unique key

    return [ key, idNormal, target, normalBam, normalBai, intervalBed ]
}.map{
    key, idNormal, target, normalBam, normalBai, intervalBed ->
    tuple (
         groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
         idNormal, target, normalBam, normalBai, intervalBed
    )
}
.transpose()
.set{ mergedChannelGermline }
}

if (runSomatic){
// --- Run Mutect2
process RunMutect2 {
  tag {idTumor + "__" + idNormal + "@" + intervalBed.baseName}

  input:
    set id, idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal), file(intervalBed) from mergedChannelSomatic 
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set id, idTumor, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forMutect2Combine

  when: "mutect2" in tools && runSomatic

  script:
  mutect2Vcf = "${idTumor}__${idNormal}_${intervalBed.baseName}.vcf.gz"
  prefix = "${mutect2Vcf}".replaceFirst(".vcf.gz", "")
  """
  gatk --java-options -Xmx8g \
    Mutect2 \
    --reference ${genomeFile} \
    --intervals ${intervalBed} \
    --input ${bamTumor} \
    --tumor-sample ${idTumor} \
    --input ${bamNormal} \
    --normal-sample ${idNormal} \
    --output ${mutect2Vcf}

  gatk --java-options -Xmx8g \
    FilterMutectCalls \
    --variant ${mutect2Vcf} \
    --stats ${prefix}.Mutect2FilteringStats.tsv \
    --output ${prefix}.filtered.vcf.gz
  """
}

//Formatting the channel to be keyed by idTumor, idNormal, and target
// group by groupKey(key, intervalBed.size())
forMutect2Combine.groupTuple().set{ forMutect2Combine }

// Combine Mutect2 VCFs, bcftools
process SomaticCombineMutect2Vcf {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${idTumor}__${idNormal}/mutect2", mode: params.publishDirMode

  input:
    set id, idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcf4Combine, mutect2CombinedVcfOutput

  when: "mutect2" in tools && runSomatic

  script:
  idTumor = id.toString().split("__")[0]
  idNormal = id.toString().split("@")[0].split("__")[1]
  target = id.toString().split("@")[1]
  outfile = "${idTumor}__${idNormal}.mutect2.vcf.gz"
  """
  bcftools concat \
    --allow-overlaps \
    ${mutect2Vcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all | \
  bcftools view \
    --samples ${idNormal},${idTumor} \
    --output-type z \
    --output-file ${outfile}

  tabix --preset vcf ${outfile}
  """
}


// --- Run Delly
Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }

process SomaticDellyCall {
  tag {idTumor + "__" + idNormal + '@' + svType}

  publishDir "${outDir}/somatic/${idTumor}__${idNormal}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  
  input:
    each svType from svTypes
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), file("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi") into dellyFilter4Combine
    set file("${idTumor}__${idNormal}_${svType}.delly.vcf.gz"), file("${idTumor}__${idNormal}_${svType}.delly.vcf.gz.tbi") into dellyOutput

  when: "delly" in tools && runSomatic

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idTumor}__${idNormal}_${svType}.bcf \
    ${bamTumor} ${bamNormal}

  echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv

  delly filter \
    --filter somatic \
    --samples samples.tsv \
    --outfile ${idTumor}__${idNormal}_${svType}.filter.bcf \
    ${idTumor}__${idNormal}_${svType}.bcf


  bcftools view --output-type z ${idTumor}__${idNormal}_${svType}.filter.bcf > ${idTumor}__${idNormal}_${svType}.delly.vcf.gz
  tabix --preset vcf ${idTumor}__${idNormal}_${svType}.delly.vcf.gz
  """
}

// --- Run Manta
process SomaticRunManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/manta", mode: params.publishDirMode, pattern: "*.manta.vcf.{gz,gz.tbi}"

  input:
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.manta.vcf.gz") into manta4Combine
    set idTumor, idNormal, target, file("${outputPrefix}.manta.vcf.gz.tbi") into mantaOutput
    set idTumor, idNormal, target, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: "manta" in tools && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  options = ""
  if (params.assayType == "exome") options = "--exome"
  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --referenceFasta ${genomeFile} \
    --normalBam ${bamNormal} \
    --tumorBam ${bamTumor} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${outputPrefix}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${outputPrefix}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${outputPrefix}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    ${outputPrefix}.manta.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    ${outputPrefix}.manta.vcf.gz.tbi
  """
}

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)

dellyFilter4Combine.groupTuple(by: [0,1,2], size: 5).combine(manta4Combine, by: [0,1,2]).set{ dellyMantaCombineChannel }

// --- Process Delly and Manta VCFs 

// Merge VCFs, Delly and Manta
process SomaticMergeDellyAndManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    set idTumor, idNormal, target, file(dellyVcfs), file(dellyVcfsIndex), file(mantaFile) from dellyMantaCombineChannel

  output:
    set val("placeHolder"), idTumor, idNormal, file("${outputPrefix}.delly.manta.vcf.gz*") into dellyMantaCombinedOutput
    set val("placeHolder"), idTumor, idNormal, file("${outputPrefix}.delly.manta.vcf.gz") into dellyMantaCombined4Aggregate
    set val("placeHolder"), idTumor, idNormal, file("${outputPrefix}.delly.manta.vcf.gz.tbi") into dellyMantaCombinedTbi4Aggregate

  when: tools.containsAll(["manta", "delly"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """ 
  bcftools view \
    --samples ${idTumor},${idNormal} \
    --output-type z \
    --output-file ${outputPrefix}.manta.swap.vcf.gz \
    ${outputPrefix}.manta.vcf.gz 
    
  tabix --preset vcf ${outputPrefix}.manta.swap.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${outputPrefix}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz ${outputPrefix}.manta.swap.vcf.gz
  
  tabix --preset vcf ${outputPrefix}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    ${outputPrefix}.delly.manta.unfiltered.vcf.gz | \
  bcftools sort \
    --output-type z \
    --output-file ${outputPrefix}.delly.manta.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.manta.vcf.gz 
  """
}


// --- Run Strelka2
bams4Strelka.combine(mantaToStrelka, by: [0, 1, 2]).set{input4Strelka}

process SomaticRunStrelka2 {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/strelka2", mode: params.publishDirMode, pattern: "*.vcf.{gz,gz.tbi}"

  input:
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal), file(mantaCSI), file(mantaCSIi) from input4Strelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])

  output:
    set idTumor, idNormal, target, file('*strelka2.vcf.gz'), file('*strelka2.vcf.gz.tbi') into strelka4Combine
    set file('*strelka2.vcf.gz'), file('*strelka2.vcf.gz.tbi') into strelkaOutput

  when: tools.containsAll(["manta", "strelka2"]) && runSomatic

  script:
  options = ""
  intervals = wgsTargets
  if (params.assayType == "exome") {
    options = "--exome"
    if (target == 'agilent') intervals = agilentTargets
    if (target == 'idt') intervals = idtTargets
  }
  outputPrefix = "${idTumor}__${idNormal}"
  outfile = "${outputPrefix}.strelka2.vcf.gz"
  """
  configureStrelkaSomaticWorkflow.py \
    ${options} \
    --reportEVSFeatures \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --indelCandidates ${mantaCSI} \
    --tumorBam ${bamTumor} \
    --normalBam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${outputPrefix}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${outputPrefix}_somatic_snvs.vcf.gz.tbi

  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz Strelka_${outputPrefix}_somatic_snvs.vcf.gz | \
  bcftools reheader \
    --samples samples.txt | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --output-type z \
    --output ${outfile}

  tabix --preset vcf ${outfile}
  """
}


mutect2CombinedVcf4Combine.combine(bamns4CombineChannel, by: [0,1,2]).combine(strelka4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }

// Combined Somatic VCFs

process SomaticCombineChannel {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutectCombinedVcf), file(mutectCombinedVcfIndex), file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal), file(strelkaVcf), file(strelkaVcfIndex) from mutectStrelkaChannel
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker, referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex
    ])
    set file(exomePoN), file(wgsPoN), file(exomePoNIndex), file(wgsPoNIndex) from Channel.value([
      referenceMap.exomePoN, referenceMap.wgsPoN,
      referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.pass.vcf") into mutationMergedVcf

  when: tools.containsAll(["manta", "strelka2", "mutect2"]) && runSomatic
  
  script:
  outputPrefix = "${idTumor}__${idNormal}"
  isecDir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWgsVcf
  if (target == "wgs") {
    pon = wgsPoN
    gnomad = gnomadWgsVcf
  }
  else {
    pon = exomePoN
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > vcf.rm.header
  echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > vcf.map.header
  echo -e "##INFO=<ID=PoN,Number=1,Type=Integer,Description=\"Count in panel of normals\">" > vcf.pon.header
  echo -e "##FORMAT=<ID=alt_count_raw,Number=1,Type=Integer,Description=\"Raw alternate allele depth\">" > vcf.ad_n.header
  echo -e "##FORMAT=<ID=alt_count_raw_fwd,Number=1,Type=Integer,Description=\"Raw alternate allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=alt_count_raw_rev,Number=1,Type=Integer,Description=\"Raw alternate allele depth on reverse strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw,Number=1,Type=Integer,Description=\"Raw reference allele depth\">" > vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw_fwd,Number=1,Type=Integer,Description=\"Raw reference allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw_rev,Number=1,Type=Integer,Description=\"Raw reference allele depth on reverse strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw,Number=1,Type=Integer,Description=\"Raw total allele depth\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw_fwd,Number=1,Type=Integer,Description=\"Raw total allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw_rev,Number=1,Type=Integer,Description=\"Raw total allele depth on reverse strand\">" >> vcf.ad_n.header

  # Get set differences of variant calls:
  # 0000: MuTect2 only
  # 0001: Strelka2 only
  # 0002: MuTect2 calls shared by Strelka2
  # 0003: Strelka2 calls shared by MuTect2

  bcftools isec \
    --output-type z \
    --prefix ${isecDir} \
    ${mutectCombinedVcf} ${strelkaVcf}

  bcftools annotate \
    --annotations ${isecDir}/0003.vcf.gz \
    --include 'FILTER!=\"PASS\"' \
    --mark-sites \"+Strelka2FILTER\" \
    -k \
    --output-type z \
    --output ${isecDir}/0003.annot.vcf.gz \
    ${isecDir}/0003.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0000.vcf.gz \
    --mark-sites +MuTect2 \
    --output-type z \
    --output ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0002.vcf.gz \
    --mark-sites \"+MuTect2;Strelka2\" \
    --output-type z \
    --output ${isecDir}/0002.tmp.vcf.gz \
    ${isecDir}/0002.vcf.gz

  tabix --preset vcf ${isecDir}/0002.tmp.vcf.gz
  tabix --preset vcf ${isecDir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${isecDir}/0003.annot.vcf.gz \
    --columns +INFO,+FORMAT,Strelka2FILTER \
    --output-type z \
    --output ${isecDir}/0002.annot.vcf.gz \
    ${isecDir}/0002.tmp.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0001.vcf.gz \
    --mark-sites +Strelka2 \
    --output-type z \
    --output ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0001.vcf.gz

  tabix --preset vcf ${isecDir}/0000.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0001.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0002.annot.vcf.gz

  # Concatenate the different sets, annotate with blacklists
  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0002.annot.vcf.gz | \
  bcftools sort | \
  bcftools annotate \
    --header-lines vcf.rm.header \
    --annotations ${repeatMasker} \
    --columns CHROM,FROM,TO,RepeatMasker | \
  bcftools annotate \
    --header-lines vcf.map.header \
    --annotations ${mapabilityBlacklist} \
    --columns CHROM,FROM,TO,EncodeDacMapability \
    --output-type z \
    --output ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.vcf.gz

  # Add gnomAD annotation
  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    --output-type z \
    --output ${idTumor}.union.gnomad.vcf.gz \
    ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.gnomad.vcf.gz

  # Add PoN annotation and flanking sequence
  bcftools annotate \
    --header-lines vcf.pon.header \
    --annotations ${pon} \
    --columns PoN:=AC_Het \
    ${idTumor}.union.gnomad.vcf.gz | \
    vt annotate_indels \
    -r ${genomeFile} \
    -o ${idTumor}.union.annot.vcf -

  # Do custom filter annotation, then filter variants
  filter-vcf.py ${idTumor}.union.annot.vcf

  mv ${idTumor}.union.annot.filter.vcf ${outputPrefix}.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}__${idNormal}.filtered.vcf \
    ${idTumor}__${idNormal}.vcf

  # Add normal read count, using all reads
  GetBaseCountsMultiSample \
    --thread ${task.cpus} \
    --maq 0 \
    --fasta ${genomeFile} \
    --bam ${idTumor}:${bamTumor} \
    --bam ${idNormal}:${bamNormal} \
    --vcf ${outputPrefix}.filtered.vcf \
    --output ${outputPrefix}.genotyped.vcf 
  
  bgzip ${outputPrefix}.filtered.vcf
  bgzip ${outputPrefix}.genotyped.vcf
  tabix --preset vcf ${outputPrefix}.filtered.vcf.gz
  tabix --preset vcf ${outputPrefix}.genotyped.vcf.gz

  bcftools annotate \
    --annotations ${outputPrefix}.genotyped.vcf.gz \
    --header-lines vcf.ad_n.header \
    --columns FORMAT/alt_count_raw:=FORMAT/AD,FORMAT/ref_count_raw:=FORMAT/RD,FORMAT/alt_count_raw_fwd:=FORMAT/ADP,FORMAT/ref_count_raw_fwd:=FORMAT/RDP,FORMAT/alt_count_raw_rev:=FORMAT/ADN,FORMAT/ref_count_raw_rev:=FORMAT/RDN,FORMAT/depth_raw:=FORMAT/DP,FORMAT/depth_raw_fwd:=FORMAT/DPP,FORMAT/depth_raw_rev:=FORMAT/DPN \
    --output-type v \
    --output ${outputPrefix}.pass.vcf \
    ${outputPrefix}.filtered.vcf.gz
  """
}

// Run vcf2maf and apply custom filters
process SomaticAnnotateMaf {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${idTumor}__${idNormal}/combined_mutations/", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

  input:
    set idTumor, idNormal, target, file(vcfMerged) from mutationMergedVcf
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFile
    file("${outputPrefix}.unfiltered.maf") into unfilteredMafOutput

  when: tools.containsAll(["manta", "strelka2", "mutect2"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}.somatic"
  mutect2InfoCols = "MBQ,MFRL,MMQ,MPOS,OCM,RPA,STR,ECNT"
  strelka2InfoCols = "RU,IC,MQ,SNVSB"
  strelka2FormatCols = "FDP,SUBDP"
  formatCols = "alt_count_raw,alt_count_raw_fwd,alt_count_raw_rev,ref_count_raw,ref_count_raw_fwd,ref_count_raw_rev,depth_raw,depth_raw_fwd,depth_raw_rev"
  formatCols = formatCols + "," + strelka2FormatCols
  if (target == "wgs") {
    infoCols = "MuTect2,Strelka2,Custom_filters,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,AC,AF,AC_nfe_seu,AF_nfe_seu,AC_afr,AF_afr,AC_nfe_onf,AF_nfe_onf,AC_amr,AF_amr,AC_eas,AF_eas,AC_nfe_nwe,AF_nfe_nwe,AC_nfe_est,AF_nfe_est,AC_nfe,AF_nfe,AC_fin,AF_fin,AC_asj,AF_asj,AC_oth,AF_oth,AC_popmax,AN_popmax,AF_popmax"
    infoCols = infoCols + "," + mutect2InfoCols + "," + strelka2InfoCols
  }
  else {
    infoCols = "MuTect2,Strelka2,Custom_filters,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax"
    infoCols = infoCols + "," + mutect2InfoCols + "," + strelka2InfoCols
  }
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 10 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --retain-fmt ${formatCols} \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0
    
  python /usr/bin/oncokb_annotator/MafAnnotator.py \
    -i ${outputPrefix}.raw.maf \
    -o ${outputPrefix}.raw.oncokb.maf

  Rscript --no-init-file /usr/bin/filter-somatic-maf.R \
    --tumor-vaf ${params.somaticVariant.tumorVaf} \
    --tumor-depth ${params.somaticVariant.tumorDepth} \
    --tumor-count ${params.somaticVariant.tumorCount} \
    --normal-depth ${params.somaticVariant.normalDepth} \
    --normal-count ${params.somaticVariant.normalCount} \
    --gnomad-allele-frequency ${params.somaticVariant.gnomadAf} \
    --normal-panel-count ${params.somaticVariant.ponCount} \
    --maf-file ${outputPrefix}.raw.oncokb.maf \
    --output-prefix ${outputPrefix}
  """
}


mafFile.into{mafFile; mafFileForMutSig}

// --- Run Mutational Signatures
process RunMutationSignatures {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.mutsig.txt") into mutSig4MetaDataParser

  when: tools.containsAll(["mutect2", "manta", "strelka2", "mutsig"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  maf2cat2.R ${outputPrefix}.somatic.maf \
  ${outputPrefix}.trinucmat.txt
  tempoSig.R --cosmic_${params.cosmic} --pvalue --nperm 10000 --seed 132 ${outputPrefix}.trinucmat.txt \
  ${outputPrefix}.mutsig.txt
  """
}



// --- Run FACETS 
process DoFacets {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "*.snp_pileup.gz"
  publishDir "${outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "${tag}_OUT.txt"
  publishDir "${outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "${outputDir}/*.{Rdata,png,out,seg,txt}"

  input:
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal) from bamFiles4DoFacets
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    file("${outfile}") into snpPileupOutput
    file("${outputDir}/*") into FacetsOutput
    set val("placeHolder"), idTumor, idNormal, file("*_OUT.txt") into FacetsOutLog4Aggregate
    set val("placeHolder"), idTumor, idNormal, file("*/*_purity.seg") into FacetsPurity4Aggregate
    set val("placeHolder"), idTumor, idNormal, file("*/*_hisens.seg") into FacetsHisens4Aggregate
    set idTumor, idNormal, target, file("${outputDir}/*purity.out") into facetsPurity4LOHHLA, facetsPurity4MetaDataParser
    set idTumor, idNormal, target, file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*hisens.cncf.txt"), val("${outputDir}") into facetsForMafAnno, facetsForMafAnnoGermline
    set idTumor, idNormal, target, file("${outputDir}/*.{Rdata,png,out,seg,txt}"), file("${idTumor}__${idNormal}.snp_pileup.gz"), val("${outputDir}") into Facets4FacetsPreview
    set val("placeHolder"), idTumor, idNormal, file("*/*.*_level.txt") into FacetsArmGeneOutput
    set val("placeHolder"), idTumor, idNormal, file("*/*.arm_level.txt") into FacetsArmLev4Aggregate
    set val("placeHolder"), idTumor, idNormal, file("*/*.gene_level.txt") into FacetsGeneLev4Aggregate
    set idTumor, idNormal, target, file("*/*.qc.txt") into FacetsQC4MetaDataParser
    set idTumor, idNormal, file("*_OUT.txt") into FacetsRunSummary

  when: "facets" in tools && runSomatic

  script:
  tag = outputFacetsSubdirectory = "${idTumor}__${idNormal}"
  outfile = tag + ".snp_pileup.gz"
  outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"
  """
  touch .Rprofile

  export SNP_PILEUP=/usr/bin/snp-pileup

  Rscript /usr/bin/facets-suite/snp-pileup-wrapper.R \
    --pseudo-snps 50 \
    --vcf-file ${facetsVcf} \
    --output-prefix ${tag} \
    --normal-bam ${bamNormal} \
    --tumor-bam ${bamTumor}

  mkdir ${outputDir}

  set +e
  i=1
  seed=\$((${params.facets.seed}-1))
  attemptNumber=0

  while [ \$i -eq 1 ]
  do
  attemptNumber=\$(( attemptNumber + 1 ))
  if [ \$attemptNumber -gt 4 ]; then 
    break
  fi
  seed=\$((seed+i))
  
  Rscript /usr/bin/facets-suite/run-facets-wrapper.R \
    --cval ${params.facets.cval} \
    --snp-window-size ${params.facets.snp_nbhd} \
    --normal-depth ${params.facets.ndepth} \
    --min-nhet ${params.facets.min_nhet} \
    --purity-cval ${params.facets.purity_cval}\
    --purity-min-nhet ${params.facets.purity_min_nhet} \
    --genome ${params.facets.genome} \
    --counts-file ${outfile} \
    --sample-id ${tag} \
    --directory ${outputDir} \
    --facets-lib-path /usr/local/lib/R/site-library \
    --seed \$seed \
    --everything \
    --legacy-output T

  i=\$?
  done
  set -e

  python3 /usr/bin/summarize_project.py \
    -p ${tag} \
    -c ${outputDir}/*cncf.txt \
    -o ${outputDir}/*out \
    -s ${outputDir}/*seg
  """
}

process DoFacetsPreviewQC {
  tag {idTumor + "__" + idNormal}
  publishDir "${outDir}/somatic/${tag}/facets/${tag}/", mode: params.publishDirMode, pattern: "${idTumor}__${idNormal}.facets_qc.txt"

  input:
    set idTumor, idNormal, target, file(facetsOutputFolderFiles), file(countsFile), facetsOutputDir from Facets4FacetsPreview
  
  output:
    set idTumor, idNormal, file("${idTumor}__${idNormal}.facets_qc.txt") into FacetsPreviewOut

  when: "facets" in tools && runSomatic

  script:
  tag = "${idTumor}__${idNormal}"
  """
  mkdir -p ${facetsOutputDir} 
  facetsFitFiles=( ${facetsOutputFolderFiles.join(" ")} )
  for i in "\${facetsFitFiles[@]}" ; do 
    cp \$i ${facetsOutputDir}/\$i
  done
  echo -e "sample_id\\tsample_path\\ttumor_id" > manifest.txt 
  echo -e "${idTumor}__${idNormal}\\t\$(pwd)\\t${idTumor}" >> manifest.txt 
  gzip manifest.txt
  mkdir -p refit_watcher/bin/ refit_watcher/refit_jobs/
  R -e "facetsPreview::generate_genomic_annotations('${idTumor}__${idNormal}', '\$(pwd)/', '/usr/bin/facets-preview/tempo_config.json')"
  cp facets_qc.txt ${idTumor}__${idNormal}.facets_qc.txt
  rm ${facetsOutputDir}/*
  """

}

FacetsRunSummary.combine(FacetsPreviewOut, by:[0,1]).into{FacetsQC4Aggregate ; FacetsQC4SomaticMultiQC} // idTumor, idNormal, summaryFiles, qcFiles
FacetsQC4Aggregate.map{ idTumor, idNormal, summaryFiles, qcFiles ->
  ["placeholder",idTumor, idNormal, summaryFiles, qcFiles]
}.set{FacetsQC4Aggregate}


// Run Polysolver
process RunPolysolver {
  tag {idNormal}
  
  input:
  set idNormal, target, file(bamNormal), file(baiNormal) from bamsNormal4Polysolver

  output:
    set val("placeHolder"), idNormal, target, file("${outputPrefix}.hla.txt") into hlaOutput, hlaOutputForLOHHLA, hlaOutputForMetaDataParser

  when: "polysolver" in tools && runSomatic
  
  script:
  outputPrefix = "${idNormal}"
  outputDir = "."
  tmpDir = "${outputDir}-nf-scratch"
  """
  cp /home/polysolver/scripts/shell_call_hla_type .
  
  sed -i "171s/TMP_DIR=.*/TMP_DIR=${tmpDir}/" shell_call_hla_type 

  bash shell_call_hla_type \
    ${bamNormal} \
    Unknown \
    1 \
    hg19 \
    STDFQ \
    0 \
    ${outputDir}

  mv winners.hla.txt ${outputPrefix}.hla.txt
  """
}

// *purity.out from FACETS, winners.hla.txt from POLYSOLVER

bamsForLOHHLA.combine(facetsPurity4LOHHLA, by: [0,1,2])
	     .combine(hlaOutputForLOHHLA, by: [1,2])
	     .set{ mergedChannelLOHHLA }

// Run LOHHLA
process RunLOHHLA {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/lohhla", mode: params.publishDirMode

  input:
    set idNormal, target, idTumor, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal), file(purityOut), placeHolder, file(winnersHla) from mergedChannelLOHHLA
    set file(hlaFasta), file(hlaDat) from Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat])

  output:
    set file("*.DNA.HLAlossPrediction_CI.txt"), file("*DNA.IntegerCPN_CI.txt"), file("*.pdf"), file("*.RData") optional true into lohhlaOutput
    set val("placeHolder"), idTumor, idNormal, file("*.DNA.HLAlossPrediction_CI.txt") into predictHLA4Aggregate
    set val("placeHolder"), idTumor, idNormal, file("*DNA.IntegerCPN_CI.txt") into intCPN4Aggregate

  when: tools.containsAll(["lohhla", "polysolver", "facets"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  cat ${winnersHla} | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
  
  PURITY=\$(grep Purity *_purity.out | grep -oP "[0-9\\.]+|NA+")
  PLOIDY=\$(grep Ploidy *_purity.out | grep -oP "[0-9\\.]+|NA+")
  cat <(echo -e "tumorPurity\ttumorPloidy") <(echo -e "${idTumor}\t\$PURITY\t\$PLOIDY") > tumor_purity_ploidy.txt

  Rscript --no-init-file /lohhla/LOHHLAscript.R \
    --patientId ${outputPrefix} \
    --normalBAMfile ${bamNormal} \
    --tumorBAMfile ${bamTumor} \
    --HLAfastaLoc ${hlaFasta} \
    --HLAexonLoc ${hlaDat} \
    --CopyNumLoc tumor_purity_ploidy.txt \
    --minCoverageFilter ${params.lohhla.minCoverageFilter} \
    --hlaPath massaged.winners.hla.txt \
    --gatkDir /picard-tools \
    --novoDir /opt/conda/bin

  if [[ -f ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt ]]
  then
    sed -i "s/^${idTumor}/${outputPrefix}/g" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt
  else
    rm -rf *.DNA.HLAlossPrediction_CI.txt
  fi

  if [[ -f ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt ]]
  then
    sed -i "s/^/${outputPrefix}\t/g" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt
    sed -i "0,/^${outputPrefix}\t/s//sample\t/" ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt
  else
    rm -rf *.DNA.IntegerCPN_CI.txt
  fi

  touch ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt
  touch ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt

  mv ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.HLAlossPrediction_CI.txt ${outputPrefix}.DNA.HLAlossPrediction_CI.txt
  mv ${outputPrefix}.${params.lohhla.minCoverageFilter}.DNA.IntegerCPN_CI.txt ${outputPrefix}.DNA.IntegerCPN_CI.txt

  if find Figures -mindepth 1 | read
  then
    mv Figures/* .
    mv ${idTumor}.minCoverage_${params.lohhla.minCoverageFilter}.HLA.pdf ${outputPrefix}.HLA.pdf
  fi
  """
}


hlaOutput.combine(mafFile, by: [1,2]).set{ input4Neoantigen }

// --- Run neoantigen prediction pipeline
process RunNeoantigen {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/neoantigen/", mode: params.publishDirMode, pattern: "*.txt"

  input:
    set idNormal, target, placeHolder, file(polysolverFile), idTumor, file(mafFile) from input4Neoantigen
    set file(neoantigenCDNA), file(neoantigenCDS) from Channel.value([
      referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS
    ])

  output:
    set val("placeHolder"), idTumor, idNormal, file("${idTumor}__${idNormal}.all_neoantigen_predictions.txt") into NetMhcStats4Aggregate
    file("${idTumor}__${idNormal}.all_neoantigen_predictions.txt") into NetMhcStatsOutput
    set idTumor, idNormal, target, file("${outputDir}/${outputPrefix}.neoantigens.maf") into mafFileForMafAnno

  when: tools.containsAll(["neoantigen", "mutect2", "manta", "strelka2"]) && runSomatic

  script:

  if (workflow.profile == "juno") {
    if(mafFile.size() > 10.MB){
      task.time = { params.maxWallTime }
    }
    else if (mafFile.size() < 5.MB){
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  outputPrefix = "${idTumor}__${idNormal}"
  outputDir = "neoantigen"
  tmpDir = "${outputDir}-tmp"
  tmpDirFullPath = "\$PWD/${tmpDir}/"  // must set full path to tmp directories for netMHC and netMHCpan to work; for some reason doesn't work with /scratch, so putting them in the process workspace
  """
  export TMPDIR=${tmpDirFullPath}
  mkdir -p ${tmpDir}
  chmod 777 ${tmpDir}

  python /usr/local/bin/neoantigen/neoantigen.py \
    --config_file /usr/local/bin/neoantigen/neoantigen-docker.config \
    --sample_id ${outputPrefix} \
    --hla_file ${polysolverFile} \
    --maf_file ${mafFile} \
    --output_dir ${outputDir}

  awk 'NR==1 {printf("%s\\t%s\\n", "sample", \$0)} NR>1 {printf("%s\\t%s\\n", "${outputPrefix}", \$0) }' neoantigen/*.all_neoantigen_predictions.txt > ${outputPrefix}.all_neoantigen_predictions.txt
  """
}


facetsForMafAnno.combine(mafFileForMafAnno, by: [0,1,2]).set{ facetsMafFileSomatic }


// --- Do FACETS MAF annotation and post processing
process SomaticFacetsAnnotation {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outputPrefix}/combined_mutations/", mode: params.publishDirMode, pattern: "*.somatic.final.maf"

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileSomatic

  output:
    set val("placeHolder"), idTumor, idNormal, file("${outputPrefix}.somatic.final.maf") into finalMaf4Aggregate
    file("file-size.txt") into mafSize
    file("${outputPrefix}.somatic.final.maf") into finalMafOutput
    set idTumor, idNormal, target, file("${outputPrefix}.somatic.final.maf") into maf4MetaDataParser

  when: tools.containsAll(["facets", "mutect2", "manta", "strelka2"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  if [ \$( cat ${maf} | wc -l ) -gt 1 ] ; then 
  Rscript --no-init-file /usr/bin/facets-suite/annotate-maf-wrapper.R \
    --facets-output ${purity_rdata} \
    --maf-file ${maf} \
    --facets-algorithm em \
    --output ${outputPrefix}.facets.maf

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-somatic.R ${outputPrefix}.facets.maf ${outputPrefix}.facets.zygosity.maf

  echo -e "${outputPrefix}\t`wc -l ${outputPrefix}.facets.zygosity.maf | cut -d ' ' -f1`" > file-size.txt

  mv ${outputPrefix}.facets.zygosity.maf ${outputPrefix}.somatic.final.maf
  else
    cp ${maf} ${outputPrefix}.somatic.final.maf
    echo -e "${outputPrefix}\t0" > file-size.txt
  fi
  """
}


// --- Run MSIsensor
process RunMsiSensor {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal)  from bamsForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.msisensor.tsv") into msi4MetaDataParser

  when: "msisensor" in tools && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  msisensor msi \
    -d ${msiSensorList} \
    -t ${bamTumor} \
    -n ${bamNormal} \
    -o ${outputPrefix}.msisensor.tsv
  """
}


facetsPurity4MetaDataParser.combine(maf4MetaDataParser, by: [0,1,2])
			   .combine(FacetsQC4MetaDataParser, by: [0,1,2])
			   .combine(msi4MetaDataParser, by: [0,1,2])
			   .combine(mutSig4MetaDataParser, by: [0,1,2])
			   .combine(hlaOutputForMetaDataParser, by: [1,2])
			   .unique()
			   .set{ mergedChannelMetaDataParser }

// --- Generate sample-level metadata
process MetaDataParser {
  tag {idTumor + "__" + idNormal}
 
  publishDir "${outDir}/somatic/${idTumor}__${idNormal}/meta_data/", mode: params.publishDirMode, pattern: "*.sample_data.txt"

  input:
    set idNormal, target, idTumor, file(purityOut), file(mafFile), file(qcOutput), file(msifile), file(mutSig), placeHolder, file(polysolverFile) from mergedChannelMetaDataParser
    set file(idtCodingBed), file(agilentCodingBed), file(wgsCodingBed) from Channel.value([
      referenceMap.idtCodingBed, referenceMap.agilentCodingBed, referenceMap.wgsCodingBed
    ]) 

  output:
    file("*.sample_data.txt") into MetaDataOutput
    set val("placeHolder"), idTumor, idNormal, file("*.sample_data.txt") into MetaData4Aggregate

  when: runSomatic

  script:
  if (target == "idt") {
    codingRegionsBed = "${idtCodingBed}"
  }
  else if (target == "agilent") {
    codingRegionsBed = "${agilentCodingBed}"
  }
  else if (target == "wgs") {
    codingRegionsBed = "${wgsCodingBed}"
  }
  """
  create_metadata_file.py \
    --sampleID ${idTumor}__${idNormal} \
    --tumorID ${idTumor} \
    --normalID ${idNormal} \
    --facetsPurity_out ${purityOut} \
    --facetsQC ${qcOutput} \
    --MSIsensor_output ${msifile} \
    --mutational_signatures_output ${mutSig} \
    --polysolver_output ${polysolverFile} \
    --MAF_input ${mafFile} \
    --coding_baits_BED ${codingRegionsBed}
  
  mv ${idTumor}__${idNormal}_metadata.txt ${idTumor}__${idNormal}.sample_data.txt
  """
}
} else { 
  if (params.pairing) {
    pairing4QC.map{ idTumor, idNormal -> 
      ["placeHolder",idTumor, idNormal,"",""]
    }.set{FacetsQC4Aggregate}
  }
}

/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/
if (runGermline){
// GATK HaplotypeCaller
process GermlineRunHaplotypecaller {
  tag {idNormal + "@" + intervalBed.baseName}

  input:
    set id, idNormal, target, file(bamNormal), file(baiNormal), file(intervalBed) from mergedChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set id, idNormal, target, file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz.tbi"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz.tbi") into haplotypecaller4Combine

  when: 'haplotypecaller' in tools && runGermline

  script:
  """
  gatk --java-options -Xmx8g \
    HaplotypeCaller \
    --reference ${genomeFile} \
    --intervals ${intervalBed} \
    --input ${bamNormal} \
    --output ${idNormal}_${intervalBed.baseName}.vcf.gz

  gatk SelectVariants \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.vcf.gz \
    --select-type-to-include SNP \
    --output ${idNormal}_${intervalBed.baseName}.snps.vcf.gz

  gatk SelectVariants \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.vcf.gz \
    --select-type-to-include INDEL \
    --output ${idNormal}_${intervalBed.baseName}.indels.vcf.gz

  gatk VariantFiltration \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.snps.vcf.gz \
    --filter-expression \"QD < 2.0\" --filter-name \"QD2\" \
    --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" \
    --filter-expression \"SOR > 3.0\" --filter-name \"SOR3\" \
    --filter-expression \"FS > 60.0\" --filter-name \"FS60\" \
    --filter-expression \"MQ < 40.0\" --filter-name \"MQ40\" \
    --filter-expression \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
    --filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
    --output ${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz

  gatk VariantFiltration \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.indels.vcf.gz \
    --filter-expression \"QD < 2.0\" --filter-name \"QD2\" \
    --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" \
    --filter-expression \"FS > 200.0\" --filter-name \"FS200\" \
    --filter-expression \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum-20\" \
    --output ${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz
  """
}

//Formatting the channel to be grouped by idTumor, idNormal, and target
// group by groupKey(key, intervalBed.size())
haplotypecaller4Combine.groupTuple().set{ haplotypecaller4Combine }

// merge VCFs, GATK HaplotypeCaller

process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  publishDir "${outDir}/germline/${idNormal}/haplotypecaller", mode: params.publishDirMode

  input:
    set id, idNormal, target, file(haplotypecallerSnpVcf), file(haplotypecallerSnpVcfIndex), file(haplotypecallerIndelVcf), file(haplotypecallerIndelVcfIndex) from haplotypecaller4Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("placeHolder"), idNormal, target, file("${outfile}"), file("${outfile}.tbi") into haplotypecallerCombinedVcf4Combine, haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools && runGermline 

  script: 
  idNormal = id.toString().split("@")[0]
  target = id.toString().split("@")[1]
  outfile = "${idNormal}.haplotypecaller.vcf.gz"  
  """
  bcftools concat \
    --allow-overlaps \
    ${haplotypecallerSnpVcf} ${haplotypecallerIndelVcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all \
    --output-type z \
    --output ${idNormal}.haplotypecaller.vcf.gz

  tabix --preset vcf ${idNormal}.haplotypecaller.vcf.gz
  """
}


// --- Run Strelka2, germline
process GermlineRunStrelka2 {
  tag {idNormal}

  publishDir "${outDir}/germline/${idNormal}/strelka2", mode: params.publishDirMode

  input:
    set idNormal, target, file(bamNormal), file(baiNormal) from bamsForStrelkaGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsIntervals),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervalsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])
    
  output:
    set val("placeHolder"), idNormal, target, file("${idNormal}.strelka2.vcf.gz"), file("${idNormal}.strelka2.vcf.gz.tbi") into strelka4CombineGermline, strelkaOutputGermline

  when: 'strelka2' in tools && runGermline
  
  script:
  options = ""
  intervals = wgsIntervals
  if (params.assayType == "exome") {
    options = "--exome"
    if (target == 'agilent') intervals = agilentTargets
    if (target == 'idt') intervals = idtTargets
  }
  """
  configureStrelkaGermlineWorkflow.py \
    ${options} \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Strelka/results/variants/variants.vcf.gz ${idNormal}.strelka2.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi ${idNormal}.strelka2.vcf.gz.tbi
  """
}

// Join HaploTypeCaller and Strelka outputs,  bcftools
haplotypecallerCombinedVcf4Combine.combine(strelka4CombineGermline, by: [0,1,2])
				  .combine(bamsTumor4VcfCombine, by: [1,2])
				  .set{ mergedChannelVcfCombine }

// --- Combine VCFs with germline calls from Haplotypecaller and Strelka2
process GermlineCombineChannel {
  tag {idTumor + "__" + idNormal}

  input:
    set idNormal, target, placeHolder, file(haplotypecallercombinedVcf), file(haplotypecallercombinedVcfIndex), file(strelkaVcf), file(strelkaVcfIndex), idTumor, file(bamTumor), file(baiTumor) from mergedChannelVcfCombine
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex,
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker, referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}__${idNormal}.germline.vcf") into mutationMergedGermline

  when: tools.containsAll(["strelka2", "haplotypecaller"]) && runGermline

  script:  
  isecDir = "${idNormal}.isec"
  gnomad = gnomadWgsVcf
  if (params.assayType == 'genome') {
    gnomad = gnomadWgsVcf
  }
  else if (params.assayType == 'exome') {
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=HaplotypeCaller,Number=0,Type=Flag,Description=\"Variant was called by HaplotypeCaller\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
  echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header

  bcftools isec \
    --output-type z \
    --prefix ${isecDir} \
    ${haplotypecallercombinedVcf} ${strelkaVcf}

  bcftools annotate \
    --annotations ${isecDir}/0003.vcf.gz \
    --include 'FILTER!=\"PASS\"' \
    --mark-sites \"+Strelka2FILTER\" \
    -k \
    --output-type z \
    --output ${isecDir}/0003.annot.vcf.gz \
    ${isecDir}/0003.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0000.vcf.gz \
    --mark-sites +HaplotypeCaller \
    --output-type z \
    --output ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0002.vcf.gz \
    --mark-sites \"+HaplotypeCaller;Strelka2\" \
    --output-type z \
    --output ${isecDir}/0002.tmp.vcf.gz \
    ${isecDir}/0002.vcf.gz

  tabix --preset vcf ${isecDir}/0002.tmp.vcf.gz
  tabix --preset vcf ${isecDir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${isecDir}/0003.annot.vcf.gz \
    --columns +FORMAT,Strelka2FILTER \
    --output-type z \
    --output ${isecDir}/0002.annot.vcf.gz \
    ${isecDir}/0002.tmp.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0001.vcf.gz \
    --mark-sites +Strelka2 \
    --output-type z \
    --output ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0001.vcf.gz

  tabix --preset vcf ${isecDir}/0000.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0001.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0002.annot.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0002.annot.vcf.gz | \
  bcftools sort | \
  bcftools annotate \
    --header-lines vcf.rm.header \
    --annotations ${repeatMasker} \
    --columns CHROM,FROM,TO,RepeatMasker | \
  bcftools annotate \
    --header-lines vcf.map.header \
    --annotations ${mapabilityBlacklist} \
    --columns CHROM,FROM,TO,EncodeDacMapability \
    --output-type z \
    --output ${idNormal}.union.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.union.pass.vcf.gz \
    ${idNormal}.union.vcf.gz

  tabix --preset vcf ${idNormal}.union.pass.vcf.gz

  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    ${idNormal}.union.pass.vcf.gz | \
  bcftools filter \
    --exclude \"${params.germlineVariant.gnomadAf}\" \
    --output-type v \
    --output ${idNormal}.union.gnomad.vcf 

  GetBaseCountsMultiSample \
    --fasta ${genomeFile} \
    --bam ${idTumor}:${bamTumor} \
    --vcf ${idNormal}.union.gnomad.vcf \
    --output ${idTumor}.genotyped.vcf

  bgzip ${idNormal}.union.gnomad.vcf
  bgzip ${idTumor}.genotyped.vcf
  tabix --preset vcf ${idNormal}.union.gnomad.vcf.gz
  tabix --preset vcf ${idTumor}.genotyped.vcf.gz

  bcftools merge \
    --output ${idTumor}__${idNormal}.germline.vcf \
    --output-type v \
    ${idNormal}.union.gnomad.vcf.gz \
    ${idTumor}.genotyped.vcf.gz
  """
}

// Run vcf2maf on combined germline VCF, apply custom filters
process GermlineAnnotateMaf {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/germline/${idNormal}/combined_mutations", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

  input:
    set idTumor, idNormal, target, file(vcfMerged) from mutationMergedGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFileGermline
    file("${outputPrefix}.unfiltered.maf") into unfilteredMafOutputGermline

  when: tools.containsAll(["strelka2", "haplotypecaller"]) && runGermline

  script:
  outputPrefix = "${idTumor}__${idNormal}.germline"
  if (target == 'wgs') {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,AC,AF,AC_nfe_seu,AF_nfe_seu,AC_afr,AF_afr,AC_nfe_onf,AF_nfe_onf,AC_amr,AF_amr,AC_eas,AF_eas,AC_nfe_nwe,AF_nfe_nwe,AC_nfe_est,AF_nfe_est,AC_nfe,AF_nfe,AC_fin,AF_fin,AC_asj,AF_asj,AC_oth,AF_oth,AC_popmax,AN_popmax,AF_popmax"
  }
  else {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax"
  }
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 4 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0

  Rscript --no-init-file /usr/bin/filter-germline-maf.R \
    --normal-depth ${params.germlineVariant.normalDepth} \
    --normal-vaf ${params.germlineVariant.normalVaf} \
    --maf-file ${outputPrefix}.raw.maf \
    --output-prefix ${outputPrefix}
  """
  }


facetsForMafAnnoGermline.combine(mafFileGermline, by: [0,1,2]).set{ facetsMafFileGermline }

process GermlineFacetsAnnotation {
  tag {idNormal}

  publishDir "${outDir}/germline/${idNormal}/combined_mutations/", mode: params.publishDirMode, pattern: "*.germline.final.maf"

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileGermline

  output:
    file("${outputPrefix}.final.maf") into mafFileOutputGermline
    set val("placeHolder"), idTumor, idNormal, file("${outputPrefix}.final.maf") into mafFile4AggregateGermline

  when: tools.containsAll(["facets", "haplotypecaller", "strelka2"]) && runGermline

  script:
  outputPrefix = "${idTumor}__${idNormal}.germline"
  """
  if [ \$( cat ${maf} | wc -l ) -gt 1 ] ; then 
  Rscript --no-init-file /usr/bin/facets-suite/annotate-maf-wrapper.R \
    --facets-output ${purity_rdata} \
    --maf-file ${maf} \
    --output ${outputPrefix}.facets.maf

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-germline.R ${outputPrefix}.facets.maf ${outputPrefix}.final.maf
  else 
    cp ${maf} ${outputPrefix}.final.maf
  fi
  """
}

// --- Call germline SVs with Delly
Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypesGermline }

process GermlineDellyCall {
  tag {idNormal + '@' + svType}

  publishDir "${outDir}/germline/${idNormal}/delly", mode: params.publishDirMode, pattern: "*delly.vcf.{gz,gz.tbi}"

  input:
    each svType from svTypesGermline
    set idNormal, target, file(bamNormal), file(baiNormal) from bamsForDellyGermline
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idNormal, target, file("${idNormal}_${svType}.delly.vcf.gz"), file("${idNormal}_${svType}.delly.vcf.gz.tbi") into dellyFilter4CombineGermline
    set file("*delly.vcf.gz"), file("*delly.vcf.gz.tbi") into dellyOutputGermline

  when: 'delly' in tools && runGermline

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idNormal}_${svType}.bcf \
    ${bamNormal}

  delly filter \
    --filter germline \
    --outfile ${idNormal}_${svType}.filter.bcf \
    ${idNormal}_${svType}.bcf

  bcftools view --output-type z ${idNormal}_${svType}.filter.bcf > ${idNormal}_${svType}.delly.vcf.gz
  tabix --preset vcf ${idNormal}_${svType}.delly.vcf.gz
  """
}

// --- Run Manta, germline
process GermlineRunManta {
  tag {idNormal}

  publishDir "${outDir}/germline/${idNormal}/manta", mode: params.publishDirMode

  input:
    set idNormal, target, file(bamNormal), file(baiNormal) from bamsForMantaGermline
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idNormal, target, file("${idNormal}.manta.vcf.gz"), file("${idNormal}.manta.vcf.gz.tbi") into manta4CombineGermline, mantaOutputGermline

  when: 'manta' in tools && runGermline

  // flag with --exome if exome
  script:
  options = ""
  if (params.assayType == "exome") options = "--exome"
  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --reference ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    ${idNormal}.manta.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    ${idNormal}.manta.vcf.gz.tbi
  """
}

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idNormal`, and `target`)


dellyFilter4CombineGermline.groupTuple(by: [0,1], size: 5).combine(manta4CombineGermline, by: [0,1]).set{ dellyMantaChannelGermline }

// --- Merge Delly and Manta VCFs 
process GermlineMergeDellyAndManta {
  tag {idNormal}

  publishDir "${outDir}/germline/${idNormal}/combined_svs/", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    set idNormal, target, file(dellyVcf), file(dellyVcfIndex), file(mantaVcf), file(mantaVcfIndex) from dellyMantaChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set val("placeHolder"), val("noTumor"), idNormal, file("${idNormal}.delly.manta.vcf.gz*") into dellyMantaCombinedOutputGermline
    set val("placeHolder"), val("noTumor"), idNormal, file("${idNormal}.delly.manta.vcf.gz") into dellyMantaCombined4AggregateGermline
    set val("placeHolder"), val("noTumor"), idNormal, file("${idNormal}.delly.manta.vcf.gz.tbi") into dellyMantaCombinedTbi4AggregateGermline

  when: tools.containsAll(["manta", "delly"]) && runGermline

  script:
  """ 
  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${idNormal}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz ${idNormal}.manta.vcf.gz

  tabix --preset vcf ${idNormal}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.delly.manta.vcf.gz \
    ${idNormal}.delly.manta.unfiltered.vcf.gz 
    
  tabix --preset vcf ${idNormal}.delly.manta.vcf.gz
  """
}
}

/*
================================================================================
=                              Quality Control                                 =
================================================================================
*/

if (runQC) {
// GATK CollectHsMetrics, WES only
bamsBQSR4CollectHsMetrics.into{bamsBQSR4CollectHsMetrics; bamsBQSRSkipCollectHsMetrics}
process QcCollectHsMetrics {
  tag {idSample}

  publishDir "${outDir}/bams/${idSample}/collecthsmetrics", mode: params.publishDirMode

  input:
    set idSample, target, file(bam), file(bai) from bamsBQSR4CollectHsMetrics
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargetsList), file(agilentTargetsList), file(idtBaitsList), file(agilentBaitsList) from Channel.value([
      referenceMap.idtTargetsList, referenceMap.agilentTargetsList,
      referenceMap.idtBaitsList, referenceMap.agilentBaitsList
    ])

  output:
    set idSample, file("${idSample}.hs_metrics.txt") into collectHsMetricsOutput, collectHsMetrics4Aggregate

  when: params.assayType == "exome" && runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  memMultiplier = params.mem_per_core ? task.cpus : 1
  javaOptions = "--java-options '-Xmx" + task.memory.toString().split(" ")[0].toInteger() * memMultiplier + "g'"

  baitIntervals = ""
  targetIntervals = ""
  if (target == 'agilent'){
    baitIntervals = "${agilentBaitsList}"
    targetIntervals = "${agilentTargetsList}"
  }
  if (target == 'idt'){
    baitIntervals = "${idtBaitsList}"
    targetIntervals = "${idtTargetsList}"
  }
  """
  gatk CollectHsMetrics \
    ${javaOptions} \
    --TMP_DIR ./ \
    --INPUT ${bam} \
    --OUTPUT ${idSample}.hs_metrics.txt \
    --REFERENCE_SEQUENCE ${genomeFile} \
    --BAIT_INTERVALS ${baitIntervals} \
    --TARGET_INTERVALS ${targetIntervals}
  """
}

if (runQC && params.assayType != "exome"){
  bamsBQSRSkipCollectHsMetrics.map{ idSample, target, bam, bai ->
    [idSample, ""]
  }.into{collectHsMetricsOutput; collectHsMetrics4Aggregate}
  
}

// Alfred, BAM QC

process QcQualimap {
  tag {idSample}
  
  publishDir "${params.outDir}/bams/${idSample}/qualimap", mode: params.publishDirMode, pattern: "*.{html,tar.gz}"
  publishDir "${params.outDir}/bams/${idSample}/qualimap", mode: params.publishDirMode, pattern: "*/*"

  input:
    set idSample, target, file(bam), file(bai) from bamsBQSR4Qualimap
    set file(idtTargets), file(agilentTargets) from Channel.value([
      file(referenceMap.idtTargets.toString().replaceAll(".gz","")), file(referenceMap.agilentTargets.toString().replaceAll(".gz",""))
    ])

  output:
    set idSample, file("${idSample}_qualimap_rawdata.tar.gz") into qualimap4MultiQC, qualimap4Aggregate
    set idSample, file("*.html"), file("css/*"), file("images_qualimapReport/*") into qualimapOutput
  
  when: runQC   

  script:
  if (params.genome == "smallGRCh37"){
    idtTargets = params.genomes["GRCh37"]."idtTargets".replaceAll(".gz","")
    agilentTargets =  params.genomes["GRCh37"]."agilentTargets".replaceAll(".gz","")
  }
  if (params.assayType == "exome"){
    if ( target == "idt"){
      gffOptions = "-gff ${idtTargets}"
    } else {
      gffOptions = "-gff ${agilentTargets}"
    }
    nr = 750
    nw = 300
  } else { 
    gffOptions = "-gd HUMAN" 
    nr = 500
    nw = 300
  }
  availMem = task.cpus * task.memory.toString().split(" ")[0].toInteger()
  // javaMem = availMem > 20 ? availMem - 4 : ( availMem > 10 ? availMem - 2 : ( availMem > 1 ? availMem - 1 : 1 ))
  javaMem = availMem > 20 ? (availMem * 0.75).round() : ( availMem > 1 ? availMem - 1 : 1 )
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }
  """
  qualimap bamqc \
  -bam ${bam} \
  ${gffOptions} \
  -outdir ${idSample} \
  -nt ${ task.cpus * 2 } \
  -nw ${nw} \
  -nr ${nr} \
  --java-mem-size=${javaMem}G

  mv ${idSample}/* . 
  tar -czf ${idSample}_qualimap_rawdata.tar.gz genome_results.txt raw_data_qualimapReport/* 
  """
}

Channel.from(true, false).set{ ignore_read_groups }
process QcAlfred {
  tag {idSample + "@" + "ignore_rg_" + ignore_rg }

  publishDir "${outDir}/bams/${idSample}/alfred", mode: params.publishDirMode

  input:
    each ignore_rg from ignore_read_groups
    set idSample, target, file(bam), file(bai) from bamsBQSR4Alfred
    file(genomeFile) from Channel.value([referenceMap.genomeFile])
    set file(idtTargets), file(agilentTargets), file(idtTargetsIndex), file(agilentTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex
    ])

  output:
    set idSample, file("${idSample}.alfred*tsv.gz") into bamsQcStats4Aggregate
    set idSample, file("${idSample}.alfred*tsv.gz"), file("${idSample}.alfred*tsv.gz.pdf") into alfredOutput

  when: runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { params.maxWallTime }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.medWallTime } : { params.minWallTime }
    }
    else {
      task.time = task.exitStatus.toString() in wallTimeExitCode ? { params.maxWallTime } : { params.medWallTime }
    }
    task.time = task.attempt < 3 ? task.time : { params.maxWallTime }
  }

  options = ""
  if (params.assayType == "exome") {
    if (target == "agilent") options = "--bed ${agilentTargets}"
    if (target == "idt") options = "--bed ${idtTargets}"
  }
  def ignore = ignore_rg ? "--ignore" : ""
  def outfile = ignore_rg ? "${idSample}.alfred.tsv.gz" : "${idSample}.alfred.per_readgroup.tsv.gz"
  """
  alfred qc ${options} \
    --reference ${genomeFile} \
    ${ignore} \
    --outfile ${outfile} \
    ${bam} && \
    Rscript --no-init-file /opt/alfred/scripts/stats.R ${outfile}
  """
}

alfredOutput
  .groupTuple(size:2, by:0)
  .join(fastPJson4sampleMultiQC, by:0)
  .join(qualimap4MultiQC, by:0)
  .join(collectHsMetricsOutput, by:0)
  .set{sampleMetrics4MultiQC}

process SampleRunMultiQC {
  tag {idSample}
  label 'multiqc_process'

  publishDir "${outDir}/bams/${idSample}/multiqc", mode: params.publishDirMode  
  
  input:
    set idSample, file(alfredRGNTsvFile), file(alfredRGYTsvFile), file(fastpJsonFile), file(qualimapFolder), file(hsmetricsFile) from sampleMetrics4MultiQC
    set file("exome_multiqc_config.yaml"), file("wgs_multiqc_config.yaml"), file("tempoLogo.png") from Channel.value([multiqcWesConfig,multiqcWgsConfig,multiqcTempoLogo])

  output:
    set idSample, file("*multiqc_report*.html"), path("*multiqc_data*.zip") into sample_multiqc_report
    set idSample, file("QC_Status.txt")

  script: 
  if (params.assayType == "exome") {
    assay = "exome"
  }
  else {
    assay = 'wgs'
  }
  """
  for i in ./*_qualimap_rawdata.tar.gz ; do 
    newFolder=\$(basename \$i | rev | cut -f 3- -d. | cut -f 3- -d_ | rev ) 
    mkdir -p qualimap/\$newFolder
    tar -xzf \$i -C qualimap/\$newFolder
  done
  
  parse_alfred.py --alfredfiles *alfred*tsv.gz 
  mkdir -p ignoreFolder 
  find . -maxdepth 1 \\( -name 'CO_ignore*mqc.yaml' -o -name 'IS_*mqc.yaml' -o -name 'GC_ignore*mqc.yaml' -o -name 'ME_aware_mqc.yaml' \\) -type f -print0 | xargs -0r mv -t ignoreFolder
  if [[ "${params.assayType}" == "exome" ]] ; then 
    find . -maxdepth 1 -name 'CM_*mqc.yaml' -type f -print0 | xargs -0r mv -t ignoreFolder
  fi

  mkdir -p fastp_original ; mv *fastp.json fastp_original
  for i in fastp_original/*fastp.json ; do 
    outname=\$(basename \$i)
    python -c "import json
  with open('\$i', 'r') as data_file:
    data = json.load(data_file)
  data_file.close()
  keys = list(data.keys())
  for element in keys:
    if element in ['read1_after_filtering','read2_after_filtering'] :
      data.pop(element, None)
  keys = list(data['summary'].keys())
  for element in keys:
    if element in ['after_filtering'] :
      data['summary'].pop(element, None)
  with open('\$outname', 'w') as data_file:
    json.dump(data, data_file)
  data_file.close()"
  done

  echo -e "\\tCoverage" > coverage_split.txt
  cover=\$(grep -i "mean cover" ./${idSample}/genome_results.txt | cut -f 2 -d"=" | sed "s/\\s*//g" | tr -d "X")
  echo -e "${idSample}\\t\${cover}" >> coverage_split.txt
  
  cp ${assay}_multiqc_config.yaml multiqc_config.yaml

  multiqc . -x ignoreFolder/ -x fastp_original/
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  multiqc . --cl_config "title: \\"Sample MultiQC Report\\"" --cl_config "subtitle: \\"${idSample} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "report_comment: \\"This report includes FASTQ and alignment statistics for the sample ${idSample}.<br/>This report does not include QC metrics from the Tumor/Normal pair that includes ${idSample}. To review pairing QC, please refer to the multiqc_report.html from the somatic-level folder.<br/>To review qc from all samples and Tumor/Normal pairs from a cohort in a single report, please refer to the multiqc_report.html from the cohort-level folder.\\"" -t "tempo" -z -x ignoreFolder/ -x fastp_original/
  mv genstats-QC_Status.txt QC_Status.txt
  """

}

//doing QcPileup and QcConpair/QcConpairAll only when --pairing [tsv] is given

if (pairingQc) {
process QcPileup {
  tag {idSample}

  publishDir "${outDir}/bams/${idSample}/pileup/", mode: params.publishDirMode
4
  input:
    set idSample, target, file(bam), file(bai) from bamsBQSR4QcPileup
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idSample, file("${idSample}.pileup") into pileupOutput, tumorPileups, normalPileups

  when: "pileup" in tools && runQC

  script:
  gatkPath = "/usr/bin/GenomeAnalysisTK.jar"
  conpairPath = "/usr/bin/conpair"
  markersBed = ""
  if (params.genome == "GRCh37") {
    markersBed = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed"
  }
  else {
    markersBed = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
  }

  if (params.mem_per_core) {
    mem = task.memory.toString().split(" ")[0].toInteger() - 1
  }
  else {
    mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
  }
  javaMem = "${mem}g"
  """
  ${conpairPath}/scripts/run_gatk_pileup_for_sample.py \
    --gatk=${gatkPath} \
    --bam=${bam} \
    --markers=${markersBed} \
    --reference=${genomeFile} \
    --xmx_java=${javaMem} \
    --outfile=${idSample}.pileup
  """
}

pairingTN.into{pairingT4Conpair; pairingN4Conpair}

tumorPileups.combine(pairingT4Conpair)
                        .filter { item ->
                          def idSample = item[0]
                          def samplePileup = item[1]
                          def idTumor = item[2]
                          def idNormal = item[3]
                          idSample == idTumor
                        }.map { item ->
                          def idTumor = item[2]
                          def idNormal = item[3]
                          def tumorPileup = item[1]
                          return [ idTumor, idNormal, tumorPileup ]
                        }
			.unique()
			.into{pileupT; pileupT4Combine}

normalPileups.combine(pairingN4Conpair)
                        .filter { item ->
                          def idSample = item[0]
                          def samplePileup = item[1]
                          def idTumor = item[2]
                          def idNormal = item[3]
                          idSample == idNormal
                        }.map { item ->
                          def idTumor = item[2]
                          def idNormal = item[3]
                          def normalPileup = item[1]
                          return [ idTumor, idNormal, normalPileup ]
                        }
			.unique()
			.into{pileupN; pileupN4Combine}


pileupT.combine(pileupN, by: [0, 1]).unique().set{ pileupConpair }

process QcConpair {
  tag {idTumor + "__" + idNormal}

  publishDir "${outDir}/somatic/${outPrefix}/conpair/", mode: params.publishDirMode

  input:
    set idTumor, idNormal, file(pileupTumor), file(pileupNormal) from pileupConpair
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.{concordance,contamination}.txt") into conpairOutput
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.concordance.txt") into conpairConcord
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.contamination.txt") into conpairContami

  when: "conpair" in tools && runQC

  script:
  outPrefix = "${idTumor}__${idNormal}"
  conpairPath = "/usr/bin/conpair"

  markersTxt = ""
  if (params.genome == "GRCh37") {
    markersTxt = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
  }
  else {
    markersTxt = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
  }

  """
  touch .Rprofile # calls to R inside the python scripts make this necessary to avoid loading user .Rprofile
  
  # Make pairing file
  echo "${idNormal}\t${idTumor}" > pairing.txt

   # Verify concordance
  ${conpairPath}/scripts/verify_concordances.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --normal_homozygous_markers_only \
    --outpre=${outPrefix}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${outPrefix}

  mv ${outPrefix}_concordance.txt ${outPrefix}.concordance.txt
  mv ${outPrefix}_contamination.txt ${outPrefix}.contamination.txt
  """
}

if (runSomatic) {
  conpairOutput
    .map{ placeHolder, idTumor, idNormal, conpairFiles -> 
      [idTumor, idNormal, conpairFiles]
    }.join(FacetsQC4SomaticMultiQC, by:[0,1])
    .set{ somaticMultiQCinput }
} else {
  conpairOutput
    .map{ placeHolder, idTumor, idNormal, conpairFiles -> 
      [idTumor, idNormal, conpairFiles, "", ""]
    }.set{ somaticMultiQCinput }

}

process SomaticRunMultiQC {
   tag {idTumor + "__" + idNormal}
   label 'multiqc_process'

  publishDir "${outDir}/somatic/${outPrefix}/multiqc", mode: params.publishDirMode  
  
  input:
    set idTumor, idNormal, file(conpairFiles), file(facetsSummaryFiles), file(facetsQCFiles) from somaticMultiQCinput
    set file("exome_multiqc_config.yaml"), file("wgs_multiqc_config.yaml"), file("tempoLogo.png") from Channel.value([multiqcWesConfig,multiqcWgsConfig,multiqcTempoLogo])

  output:
    set idTumor, idNormal, file("*multiqc_report*.html"), file("*multiqc_data*.zip") into somatic_multiqc_report

  script: 
  outPrefix = "${idTumor}__${idNormal}"
  if (params.assayType == "exome") {
    assay = "exome"
  }
  else {
    assay = 'wgs'
  }
  """
  echo -e "\\tTumor\\tNormal\\tTumor_Contamination\\tNormal_Contamination\\tConcordance" > conpair.tsv
  for i in ./*contamination.txt ; do 
     j=./\$(basename \$i | cut -f 1 -d.).concordance.txt
     echo -e "\$(tail -n +2 \$i | sort -r | cut -f 2| head -1)\\t\$(tail -n +2 \$i | sort -r | cut -f 2| paste -sd"\\t")\\t\$(tail -n +2 \$i | sort -r | cut -f 3| paste -sd"\\t")\\t\$(tail -1 \$j | cut -f 2 )" >> conpair.tsv
  done

  head -1 ${facetsQCFiles} | cut -f 1,28,97  | sed "s/^tumor_sample_id//g"> ${facetsQCFiles}.qc.txt
  tail -n +2 ${facetsQCFiles} | cut -f 1,28,97 | sed "s/TRUE\$/PASS/g" | sed "s/FALSE\$/FAIL/g" >> ${facetsQCFiles}.qc.txt

  cp conpair.tsv conpair_genstat.tsv
  mkdir -p ignoreFolder ; mv conpair.tsv ignoreFolder
  cp ${assay}_multiqc_config.yaml multiqc_config.yaml
  
  multiqc . -x ignoreFolder
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  multiqc . --cl_config "title: \\"Somatic MultiQC Report\\"" --cl_config "subtitle: \\"${outPrefix} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "report_comment: \\"This report includes QC statistics related to the Tumor/Normal pair ${outPrefix}.<br/>This report does not include FASTQ or alignment QC of either ${idTumor} or ${idNormal}. To review FASTQ and alignment QC, please refer to the multiqc_report.html from the bam-level folder.<br/>To review qc from all samples and Tumor/Normal pairs from a cohort in a single report, please refer to the multiqc_report.html from the cohort-level folder.\\"" -z -x ignoreFolder

  """

}

if(runConpairAll){
pileupT4Combine.combine(pileupN4Combine).unique().set{ pileupConpairAll }

process QcConpairAll {
  tag {idTumor + "@" + idNormal}

  input:
    set idTumor, idNormal_noUse, file(pileupTumor), idTumor_noUse, idNormal, file(pileupNormal) from pileupConpairAll
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.{concordance,contamination}.txt") into conpairAllOutput
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.concordance.txt") into conpairAllConcord
    set val("placeHolder"), idTumor, idNormal, file("${outPrefix}.contamination.txt") into conpairAllContami

  when: runConpairAll && runQC

  script:
  outPrefix = "${idTumor}__${idNormal}"
  conpairPath = "/usr/bin/conpair"

  markersTxt = ""
  if (params.genome == "GRCh37") {
    markersTxt = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
  }
  else {
    markersTxt = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
  }

  """
  touch .Rprofile # calls to R inside the python scripts make this necessary to avoid loading user .Rprofile
  
  # Make pairing file
  echo "${idNormal}\t${idTumor}" > pairing.txt

   # Verify concordance
  ${conpairPath}/scripts/verify_concordances.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --normal_homozygous_markers_only \
    --outpre=${outPrefix}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${outPrefix}

  mv ${outPrefix}_concordance.txt ${outPrefix}.concordance.txt
  mv ${outPrefix}_contamination.txt ${outPrefix}.contamination.txt
  """
}
}

// -- Run based on QcConpairAll channels or the single QcConpair channels
conpairConcord4Aggregate = (!runConpairAll ? conpairConcord : conpairAllConcord)
conpairContami4Aggregate = (!runConpairAll ? conpairContami : conpairAllContami)
} // End of "if (pairingQc){}", doing QcPileup or QcConpair/QcConpairAll only when --pairing [tsv] is given

} // End of "if (runQc){}"

/*
================================================================================
=                              Cohort Aggregation                              =
================================================================================
*/

if ( !params.mapping && !params.bamMapping ){
  runSomatic = true
  runGermline = true
  runQC = true
  pairingQc = true
  if(!params.watch){
    TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
	      .set{inputAggregate}
  }
  else{
    watchAggregateWithPath(file(runAggregate, checkIfExists: true))
	   .set{inputAggregate}
    epochMap[runAggregate] = 0
  }
  inputAggregate.multiMap{ cohort, idTumor, idNormal, path ->
		  finalMaf4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.final.maf" )]
                  NetMhcStats4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.all_neoantigen_predictions.txt")]
                  FacetsPurity4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*_purity.seg")]
                  FacetsHisens4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*_hisens.seg")]
                  FacetsOutLog4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*_OUT.txt")]
                  FacetsArmLev4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*.arm_level.txt")]
                  FacetsGeneLev4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*.gene_level.txt")]
                  FacetsQC4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*_OUT.txt"), file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*.facets_qc.txt")]
                  dellyMantaCombined4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.delly.manta.vcf.gz")]
                  dellyMantaCombinedTbi4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.delly.manta.vcf.gz.tbi")]
                  predictHLA4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.DNA.HLAlossPrediction_CI.txt")]
                  intCPN4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*DNA.IntegerCPN_CI.txt")]
                  MetaData4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.sample_data.txt")]
                  mafFile4AggregateGermline: [idTumor, idNormal, cohort, "placeHolder", file(path + "/germline/" + idNormal + "/*/" + idTumor + "__" + idNormal + ".germline.final.maf")]
                  dellyMantaCombined4AggregateGermline: [idNormal, cohort, idTumor, "placeHolder", "noTumor", file(path + "/germline/" + idNormal + "/*/*.delly.manta.vcf.gz")]
                  dellyMantaCombinedTbi4AggregateGermline: [idNormal, cohort, idTumor, "placeHolder", "noTumor", file(path + "/germline/" + idNormal + "/*/*.delly.manta.vcf.gz.tbi")]
                  conpairConcord4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/" + idTumor + "__" + idNormal + ".concordance.txt")]
                  conpairContami4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/" + idTumor + "__" + idNormal + ".contamination.txt")]
		  alfredIgnoreYTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.alfred.tsv.gz/")]
		  alfredIgnoreYNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.alfred.tsv.gz/")]
		  alfredIgnoreNTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.alfred.per_readgroup.tsv.gz/")]
		  alfredIgnoreNNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.alfred.per_readgroup.tsv.gz/")]
      qualimapTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/qualimap/${idTumor}_qualimap_rawdata.tar.gz")]
      qualimapNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/qualimap/${idNormal}_qualimap_rawdata.tar.gz")]
		  hsMetricsTumor: params.assayType == "exome" ? [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.hs_metrics.txt")] : [cohort, idTumor, idNormal, ""]
		  hsMetricsNormal: params.assayType == "exome" ? [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.hs_metrics.txt")] : [cohort, idTumor, idNormal, ""]
      fastpTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.fastp.json")]
      fastpNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.fastp.json")]
                }
		.set {aggregateList}
  inputSomaticAggregateMaf = aggregateList.finalMaf4Aggregate.transpose().groupTuple(by:[2])
  inputSomaticAggregateNetMHC = aggregateList.NetMhcStats4Aggregate.transpose().groupTuple(by:[2])
  inputPurity4Aggregate = aggregateList.FacetsPurity4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputHisens4Aggregate = aggregateList.FacetsHisens4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputOutLog4Aggregate = aggregateList.FacetsOutLog4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputArmLev4Aggregate = aggregateList.FacetsArmLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputGeneLev4Aggregate = aggregateList.FacetsGeneLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputFacetsQC4CohortMultiQC = aggregateList.FacetsQC4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
  inputSomaticAggregateSv = aggregateList.dellyMantaCombined4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputSomaticAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputPredictHLA4Aggregate = aggregateList.predictHLA4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputIntCPN4Aggregate = aggregateList.intCPN4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
  inputSomaticAggregateMetadata = aggregateList.MetaData4Aggregate.transpose().groupTuple(by:[2])
  inputGermlineAggregateMaf = aggregateList.mafFile4AggregateGermline.transpose().groupTuple(by:[2])
  inputGermlineAggregateSv = aggregateList.dellyMantaCombined4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
  inputGermlineAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
  inputAlfredIgnoreY = aggregateList.alfredIgnoreYTumor.unique().combine(aggregateList.alfredIgnoreYNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
  inputAlfredIgnoreN = aggregateList.alfredIgnoreNTumor.unique().combine(aggregateList.alfredIgnoreNNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
  inputQualimap4CohortMultiQC = aggregateList.qualimapTumor.unique().combine(aggregateList.qualimapNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
  inputHsMetrics = aggregateList.hsMetricsTumor.unique().combine(aggregateList.hsMetricsNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
  aggregateList.conpairConcord4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.into{inputConpairConcord4Aggregate; inputConpairConcord4MultiQC}
  aggregateList.conpairContami4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.into{inputConpairContami4Aggregate; inputConpairContami4MultiQC}
  aggregateList.fastpTumor.unique().combine(aggregateList.fastpNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}.set{inputFastP4MultiQC}

}
else if(!(runAggregate == false)) {
  if (!(runAggregate == true)){
    if(!params.watch){
      TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
	        .groupTuple()
                .map{ cohort, idTumor, idNormal, pathNoUse
                      -> tuple( groupKey(cohort, idTumor instanceof Collection ? idTumor.size() : 1), idTumor, idNormal)
                }
                .transpose()
	        .set{inputAggregate}
    }
    else{
      watchAggregate(file(runAggregate, checkIfExists: false))
	     .set{inputAggregate}
      epochMap[runAggregate] = 0
    }
  }
  else if(runAggregate == true){
    inputPairing.into{inputPairing; cohortTable}
    cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
	       .set{inputAggregate}
  }
  else{}
  inputAggregate.into{ cohortSomaticAggregateMaf;
                       cohortSomaticAggregateNetMHC;
                       cohortSomaticAggregateFacets;
                       cohortSomaticAggregateFacets1;
                       cohortSomaticAggregateFacets2;
                       cohortSomaticAggregateFacets3;
                       cohortSomaticAggregateFacets4;
                       cohortSomaticAggregateSv;
                       cohortSomaticAggregateSv1;
                       cohortSomaticAggregateLOHHLA;
                       cohortSomaticAggregateLOHHLA1;
                       cohortSomaticAggregateMetadata;
                       cohortGermlineAggregateMaf;
                       cohortGermlineAggregateSv;
                       cohortGermlineAggregateSv1;
                       cohortQcBamAggregate;
                       cohortQcBamAggregate1;
                       cohortQcBamAggregate2;
                       cohortQcBamAggregate3;
                       cohortQcConpairAggregate;
                       cohortQcConpairAggregate1;
                       cohortQcFastPAggregate;
                       cohortQcFacetsAggregate
                }
  if (runSomatic){
  inputSomaticAggregateMaf = cohortSomaticAggregateMaf.combine(finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
  inputSomaticAggregateNetMHC = cohortSomaticAggregateNetMHC.combine(NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
  inputPurity4Aggregate = cohortSomaticAggregateFacets.combine(FacetsPurity4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputHisens4Aggregate = cohortSomaticAggregateFacets1.combine(FacetsHisens4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputOutLog4Aggregate = cohortSomaticAggregateFacets2.combine(FacetsOutLog4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputArmLev4Aggregate = cohortSomaticAggregateFacets3.combine(FacetsArmLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputGeneLev4Aggregate = cohortSomaticAggregateFacets4.combine(FacetsGeneLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputSomaticAggregateSv = cohortSomaticAggregateSv.combine(dellyMantaCombined4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputSomaticAggregateSvTbi = cohortSomaticAggregateSv1.combine(dellyMantaCombinedTbi4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputPredictHLA4Aggregate = cohortSomaticAggregateLOHHLA.combine(predictHLA4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputIntCPN4Aggregate = cohortSomaticAggregateLOHHLA1.combine(intCPN4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
  inputSomaticAggregateMetadata = cohortSomaticAggregateMetadata.combine(MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])

  if (runGermline){
  inputGermlineAggregateMaf = cohortGermlineAggregateMaf.combine(mafFile4AggregateGermline, by:[1,2]).groupTuple(by:[2])
  inputGermlineAggregateSv = cohortGermlineAggregateSv.combine(dellyMantaCombined4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
  inputGermlineAggregateSvTbi = cohortGermlineAggregateSv1.combine(dellyMantaCombinedTbi4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
  }
  }

  if (runQC){
  inputPairing.into{inputPairing;inputPairing1;inputPairing2;inputPairing3;inputPairing4;inputPairing5}
  bamsQcStats4Aggregate.branch{ item ->
  			def idSample = item[0]
  			def alfred = item[1]
  
  			ignoreY: alfred =~ /.+\.alfred\.tsv\.gz/
  			ignoreN: alfred =~ /.+\.alfred\.per_readgroup\.tsv\.gz/
  		     }
  		     .set{bamsQcStats4Aggregate}

  inputPairing1.combine(bamsQcStats4Aggregate.ignoreY)
  	     .branch { item ->
  		def idTumor = item[0]
  		def idNormal = item[1]
  		def idSample = item[2]
  		def alfred = item[3]
  
  		tumor: idSample == idTumor
  		normal: idSample == idNormal
  	    }
  	    .set{alfredIgnoreY}
  alfredIgnoreY.tumor.combine(alfredIgnoreY.normal, by:[0,1])
  		   .combine(cohortQcBamAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
  		   .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
  		   .groupTuple(by:[0])
  		   .map{ item ->
  			def cohort = item[0]
  			def idTumors = item[1].unique()
  			def idNormals = item[2].unique()
  			def fileTumor = item[3].unique()
  			def fileNormal = item[4].unique()
  
			  [cohort, fileTumor, fileNormal]
  		   }
  		   .unique()
  		   .set{alfredIgnoreY}
  
  inputPairing2.combine(bamsQcStats4Aggregate.ignoreN)
  	     .branch { item ->
  		def idTumor = item[0]
  		def idNormal = item[1]
  		def idSample = item[2]
  		def alfred = item[3]
  
  		tumor: idSample == idTumor
  		normal: idSample == idNormal
  	    }
  	    .set{alfredIgnoreN}
  alfredIgnoreN.tumor.combine(alfredIgnoreN.normal, by:[0,1])
  		   .combine(cohortQcBamAggregate1.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
  		   .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
  		   .groupTuple(by:[0])
  		   .map{ item ->
  			def cohort = item[0]
  			def idTumors = item[1].unique()
  			def idNormals = item[2].unique()
  			def fileTumor = item[3].unique()
  			def fileNormal = item[4].unique()
  
			[cohort, fileTumor, fileNormal]
  		   }
  		   .unique()
  		   .set{alfredIgnoreN}
  
  inputPairing3.combine(collectHsMetrics4Aggregate)
  	     .branch { item ->
  		def idTumor = item[0]
  		def idNormal = item[1]
  		def idSample = item[2]
  		def hsMetrics = item[3]
  
  		tumor: idSample == idTumor
  		normal: idSample == idNormal
  	    }
  	    .set{hsMetrics}
  hsMetrics.tumor.combine(hsMetrics.normal, by:[0,1])
  	       .combine(cohortQcBamAggregate2.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
  	       .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
  	       .groupTuple(by:[0])
  	       .map{ item ->
  		    def cohort = item[0]
  		    def idTumors = item[1].unique()
  		    def idNormals = item[2].unique()
  		    def fileTumor = item[3].unique()
  		    def fileNormal = item[4].unique()
  		    
		    [cohort, fileTumor, fileNormal]
  	       }
  	       .unique()
  	       .set{hsMetrics}
  inputHsMetrics = hsMetrics
  inputPairing4.combine(fastPJson4cohortMultiQC)
         .branch { item ->
      def idTumor = item[0]
      def idNormal = item[1]
      def idSample = item[2]
      def jsonFiles = item[3]
  
      tumor: idSample == idTumor
      normal: idSample == idNormal
        }
        .set{fastPMetrics}
  fastPMetrics.tumor.combine(fastPMetrics.normal, by:[0,1])
           .combine(cohortQcFastPAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
           .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
           .groupTuple(by:[0])
           .map{ item ->
          def cohort = item[0]
          def idTumors = item[1].unique()
          def idNormals = item[2].unique()
          def fileTumor = item[3].flatten().unique()
          def fileNormal = item[4].flatten().unique()
          
          [cohort, fileTumor, fileNormal]
           }
           .unique()
           .set{fastPMetrics}

  inputPairing5.combine(qualimap4Aggregate)
    .branch { idTumor, idNormal, idSample, qualimapDir ->
      tumor:  idSample == idTumor
      normal: idSample == idNormal
    }
    .set{qualimap4AggregateTN}
  qualimap4AggregateTN.tumor.combine(qualimap4AggregateTN.normal, by:[0,1])
    .combine(cohortQcBamAggregate3.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
    .map{ item -> [item[6], item[3], item[5]]}
    .groupTuple(by:[0])
    .map{ cohort, fileTumor, fileNormal ->
        [cohort, fileTumor.unique(), fileNormal.unique()]
    }
    .unique()
    .set{inputQualimap4CohortMultiQC}

  inputAlfredIgnoreY = alfredIgnoreY
  inputAlfredIgnoreN = alfredIgnoreN
  inputFastP4MultiQC = fastPMetrics
  if (pairingQc){
  cohortQcConpairAggregate.combine(conpairConcord4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.into{inputConpairConcord4Aggregate; inputConpairConcord4MultiQC}
  cohortQcConpairAggregate1.combine(conpairContami4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.into{inputConpairContami4Aggregate; inputConpairContami4MultiQC}
  cohortQcFacetsAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}.set{inputFacetsQC4CohortMultiQC}
  }
  }
}
else{}

if (runAggregate && runSomatic) {
process SomaticAggregateMaf {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set idTumors, idNormals, cohort, placeHolder, file(mafFile) from inputSomaticAggregateMaf

  output:
    file("mut_somatic.maf") into mutationAggregatedOutput

  when: runSomatic

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge MAF files
  mkdir mut
  mv *.maf mut/
  for i in mut/*.maf ; do 
    if [ \$( cat \$i | wc -l ) -gt 1 ] ; then 
      cat \$i
    fi
  done | grep ^Hugo_Symbol | head -n 1 > mut_somatic.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo_Symbol" | sort -k5,5V -k6,6n >> mut_somatic.maf
  """
}

process SomaticAggregateNetMHC {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set idTumors, idNormals, cohort, placeHolder, file(netmhcCombinedFile) from inputSomaticAggregateNetMHC

  output:
    file("mut_somatic_neoantigens.txt") into NetMhcAggregatedOutput

  when: runSomatic

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp
  ## Collect and merge neoantigen prediction
  mkdir neoantigen
  mv *.all_neoantigen_predictions.txt neoantigen/
  awk 'FNR==1 && NR!=1{next;}{print}' neoantigen/*.all_neoantigen_predictions.txt > mut_somatic_neoantigens.txt
  """
}

inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
		     .join(inputOutLog4Aggregate, by:[0])
		     .join(inputArmLev4Aggregate, by:[0])
		     .join(inputGeneLev4Aggregate, by:[0])
		     .set{inputSomaticAggregateFacets}

process SomaticAggregateFacets {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(purity), file(Hisens), file(outLog), file(armLev), file(geneLev) from inputSomaticAggregateFacets

  output:
    set file("cna_hisens_run_segmentation.seg"), file("cna_purity_run_segmentation.seg"), file("cna_armlevel.txt"), file("cna_genelevel.txt"), file("cna_facets_run_info.txt") into FacetsAnnotationAggregatedOutput

  when: runSomatic

  script:
  """
  # Collect and merge FACETS outputs
  # Arm-level and gene-level output is filtered
  mkdir facets_tmp
  mv *_OUT.txt facets_tmp/
  mv *{purity,hisens}.seg facets_tmp/

  awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_hisens.seg > cna_hisens_run_segmentation.seg
  awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_purity.seg > cna_purity_run_segmentation.seg
  awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_OUT.txt > cna_facets_run_info.txt
  mv *{gene_level,arm_level}.txt facets_tmp/
  cat facets_tmp/*gene_level.txt | head -n 1 > cna_genelevel.txt
  awk -v FS='\t' '{ if (\$24 != "DIPLOID" && (\$25 == "PASS" || \$25 == "RESCUE" ))  print \$0 }' facets_tmp/*gene_level.txt >> cna_genelevel.txt
  cat facets_tmp/*arm_level.txt | head -n 1 > cna_armlevel.txt
  cat facets_tmp/*arm_level.txt | grep -v "DIPLOID" | grep -v "Tumor_Sample_Barcode" >> cna_armlevel.txt || [[ \$? == 1 ]]
  """
}

inputSomaticAggregateSv.join(inputSomaticAggregateSvTbi)
		       .set{inputSomaticAggregateSv}

process SomaticAggregateSv {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(dellyMantaVcf), file(dellyMantaVcfTbi) from inputSomaticAggregateSv

  output:
    file("sv_somatic.vcf.{gz,gz.tbi}") into svAggregatedOutput

  when: runSomatic

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge Delly and Manta VCFs
  mkdir sv/
  mv *delly.manta.vcf.gz* sv/
  vcfs=(\$(ls sv/*delly.manta.vcf.gz))
  if [[ \${#vcfs[@]} > 1 ]]
  then
    bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output sv_somatic.vcf.gz \
    sv/*delly.manta.vcf.gz
  else
    mv \${vcfs[0]} sv_somatic.vcf.gz
  fi

  tabix --preset vcf sv_somatic.vcf.gz
  """
}

inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
			 .set{inputSomaticAggregateLOHHLA}

process SomaticAggregateLOHHLA {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(preditHLA), file(intCPN) from inputSomaticAggregateLOHHLA

  output:
    file("DNA.IntegerCPN_CI.txt") into lohhlaDNAIntegerCPNOutput
    file("HLAlossPrediction_CI.txt") into lohhlaHLAlossPredictionAggregatedOutput

  when: runSomatic

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp
  mkdir lohhla
  mv *.txt lohhla/
  awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*HLAlossPrediction_CI.txt > HLAlossPrediction_CI.txt
  awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*DNA.IntegerCPN_CI.txt > DNA.IntegerCPN_CI.txt
  """
}

process SomaticAggregateMetadata {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set idTumors, idNormals, cohort, placeHolder, file(metaDataFile) from inputSomaticAggregateMetadata

  output:
    file("sample_data.txt") into MetaDataAggregatedOutput

  when: runSomatic

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge metadata file
  mkdir sample_data_tmp
  mv *.sample_data.txt sample_data_tmp/
  awk 'FNR==1 && NR!=1{next;}{print}' sample_data_tmp/*.sample_data.txt > sample_data.txt
  """
}
}


if (runAggregate && runGermline) {
// --- Aggregate per-sample germline data, MAF
process GermlineAggregateMaf {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set idTumors, idNormals, cohort, placeHolder, file(mafFile) from inputGermlineAggregateMaf

  output:
    file("mut_germline.maf") into mutationAggregatedGermlineOutput

  when: runGermline

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge MAF files
  mkdir mut
  mv *.maf mut/
  for i in mut/*.maf ; do 
    if [ \$( cat \$i | wc -l ) -gt 1 ] ; then 
      cat \$i
    fi
  done | grep ^Hugo | head -n1 > mut_germline.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo" | sort -k5,5V -k6,6n >> mut_germline.maf

  """
}

// --- Aggregate per-sample germline data, SVs
inputGermlineAggregateSv.join(inputGermlineAggregateSvTbi)
			.set{inputGermlineAggregateSv}

process GermlineAggregateSv {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(dellyMantaVcf), file(dellyMantaVcfTbi) from inputGermlineAggregateSv

  output:
    file("sv_germline.vcf.{gz,gz.tbi}") into svAggregatedGermlineOutput

  when: runGermline

  script:
  """
  ## Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  ## Collect and merge Delly and Manta VCFs
  mkdir sv
  mv  *.delly.manta.vcf.gz* sv/
  vcfs=(\$(ls sv/*delly.manta.vcf.gz))
  if [[ \${#vcfs[@]} > 1 ]]
  then
    bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output sv_germline.vcf.gz \
    sv/*delly.manta.vcf.gz
  else
    mv \${vcfs[0]} sv_germline.vcf.gz
  fi

  tabix --preset vcf sv_germline.vcf.gz
  """
}
}


if (runAggregate && runQC) {

inputAlfredIgnoreY.into{inputAlfredIgnoreY; inputAlfredIgnoreY4MultiQC }
inputAlfredIgnoreN.into{inputAlfredIgnoreN; inputAlfredIgnoreN4MultiQC }

inputHsMetrics.into{inputHsMetrics; inputHsMetrics4MultiQC}
inputAlfredIgnoreY.join(inputAlfredIgnoreN)
		  .join(inputHsMetrics)
		  .set{inputQcBamAggregate}

process QcBamAggregate {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(alfredIgnoreYTumor), file(alfredIgnoreYNormal), file(alfredIgnoreNTumor), file(alfredIgnoreNNormal), file(hsMetricsTumor), file(hsMetricsNormal) from inputQcBamAggregate

  output:
    file('alignment_qc.txt') into alignmentQcAggregatedOutput

  when: runQC

  script:
  if (params.assayType == "exome") {
    assayType = "wes"
  }
  else {
    assayType = 'wgs'
  }
  """
  Rscript --no-init-file /usr/bin/create-aggregate-qc-file.R -n ${task.cpus} -a ${assayType}
  """
}

if (pairingQc){

inputConpairConcord4Aggregate.join(inputConpairContami4Aggregate)
			     .set{inputQcConpairAggregate}

process QcConpairAggregate {

  tag {cohort}

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(concordFile), file(contamiFile) from inputQcConpairAggregate

  output:
    set file('concordance_qc.txt'), file('contamination_qc.txt') into conpairAggregatedOutput

  when: runQC

  script:
  """
  if ls *.concordance.txt 1> /dev/null 2>&1; then
    echo -e "Pair\tConcordance" > concordance_qc.txt
    grep -v "concordance" *.concordance.txt | sed 's/.concordance.txt:/\t/' | cut -f1,3 | sort -k1,1 >> concordance_qc.txt
  fi
  if ls *.contamination.txt 1> /dev/null 2>&1; then
    echo -e "Pair\tSample_Type\tSample_ID\tContamination" > contamination_qc.txt
    grep -v "Contamination" *.contamination.txt | sed 's/.contamination.txt:/\t/' | sort -k1,1 >> contamination_qc.txt
  fi
  touch concordance_qc.txt contamination_qc.txt
  """
}

inputFastP4MultiQC
  .join(inputAlfredIgnoreY4MultiQC,by:0)
  .join(inputAlfredIgnoreN4MultiQC,by:0)
  .join(inputConpairConcord4MultiQC,by:0)
  .join(inputConpairContami4MultiQC,by:0)
  .join(inputFacetsQC4CohortMultiQC,by:0)
  .join(inputQualimap4CohortMultiQC,by:0)
  .join(inputHsMetrics4MultiQC, by:0)
  .set{inputCohortRunMultiQC}

process CohortRunMultiQC {
  tag {cohort}
  label 'multiqc_process'

  publishDir "${outDir}/cohort_level/${cohort}", mode: params.publishDirMode

  input:
    set cohort, file(fastPTumor), file(fastPNormal), file(alfredIgnoreYTumor), file(alfredIgnoreYNormal), file(alfredIgnoreNTumor), file(alfredIgnoreNNormal), file(concordFile), file(contamiFile), file(FacetsSummaryFile), file(FacetsQCFile), file(qualimapFolderTumor), file(qualimapFolderNormal), file(hsMetricsTumor), file(hsMetricsNormal) from inputCohortRunMultiQC
    set file("exome_multiqc_config.yaml"), file("wgs_multiqc_config.yaml"), file("tempoLogo.png") from Channel.value([multiqcWesConfig,multiqcWgsConfig,multiqcTempoLogo])

  output:
    set cohort, file("*multiqc_report*.html"), file("*multiqc_data*.zip") into cohort_multiqc_report

  when: runQC 

  script:
  if (params.assayType == "exome") {
    assay = "exome"
  }
  else {
    assay = 'wgs'
  }
  """
  for i in ./*_qualimap_rawdata.tar.gz ; do 
    newFolder=\$(basename \$i | rev | cut -f 3- -d. | cut -f 3- -d_ | rev ) 
    mkdir -p qualimap/\$newFolder
    tar -xzf \$i -C qualimap/\$newFolder
  done
  echo -e "\\tTumor\\tNormal\\tTumor_Contamination\\tNormal_Contamination\\tConcordance" > conpair.tsv
  for i in ./*contamination.txt ; do 
     j=./\$(basename \$i | cut -f 1 -d.).concordance.txt
     echo -e "\$(tail -n +2 \$i | sort -r | cut -f 2| head -1)\\t\$(tail -n +2 \$i | sort -r | cut -f 2| paste -sd"\\t")\\t\$(tail -n +2 \$i | sort -r | cut -f 3| paste -sd"\\t")\\t\$(tail -1 \$j | cut -f 2 )" >> conpair.tsv
  done
  cp conpair.tsv conpair_genstat.tsv

  mkdir -p fastp_original ; mv *fastp.json fastp_original
  for i in fastp_original/*fastp.json ; do 
    outname=\$(basename \$i)
    python -c "import json
  with open('\$i', 'r') as data_file:
    data = json.load(data_file)
  data_file.close()
  keys = list(data.keys())
  for element in keys:
    if element in ['read1_after_filtering','read2_after_filtering'] :
      data.pop(element, None)
  keys = list(data['summary'].keys())
  for element in keys:
    if element in ['after_filtering'] :
      data['summary'].pop(element, None)
  with open('\$outname', 'w') as data_file:
    json.dump(data, data_file)
  data_file.close()"
  done
  
  for i in ./*/genome_results.txt ; do
    sampleName=\$(dirname \$i | xargs -n 1 basename )
    cover=\$(grep -i "mean cover" \$i | cut -f 2 -d"=" | sed "s/\\s*//g" | tr -d "X")
    echo -e "\${sampleName}\\t\${cover}"
  done > flatCoverage
  echo -e "\\tTumor_Coverage\\tNormal_Coverage" > coverage_split.txt
  join -1 2 -2 1 -o 1.1,1.2,1.3,2.2 -t \$'\\t' <(join -1 1 -2 1 -t \$'\\t' <(cut -f2,3 conpair_genstat.tsv | tail -n +2 | sort | uniq) <(cat flatCoverage | sort | uniq)) <(cat flatCoverage | sort | uniq) | cut -f 1,3,4 >> coverage_split.txt
  join -1 1 -2 1 -t \$'\\t' <(cut -f 3 conpair_genstat.tsv | sort | uniq | sed "s/\$/\\t/g" ) <(cat flatCoverage | sort | uniq) >> coverage_split.txt
  
  parse_alfred.py --alfredfiles *alfred*tsv.gz
  mkdir -p ignoreFolder 
  find . -maxdepth 1 \\( -name 'CO_ignore_mqc.yaml' -o -name 'IS_*mqc.yaml' -o -name 'GC_ignore_mqc.yaml' -o -name 'ME_aware_mqc.yaml' \\) -type f -print0 | xargs -0r mv -t ignoreFolder
  if [[ "${params.assayType}" == "exome" ]] ; then 
    find . -maxdepth 1 -name 'CM_*mqc.yaml' -type f -print0 | xargs -0r mv -t ignoreFolder
  fi
  mv conpair.tsv ignoreFolder

  for i in *.facets_qc.txt ; do 
    head -1 \$i | cut -f 1,28,97  | sed "s/^tumor_sample_id//g"> \$i.qc.txt
    tail -n +2 \$i | cut -f 1,28,97 | sed "s/TRUE\$/PASS/g" | sed "s/FALSE\$/FAIL/g" >> \$i.qc.txt
  done

  cp ${assay}_multiqc_config.yaml multiqc_config.yaml

  samplesNum=`for i in ./*contamination.txt ; do tail -n +2 \$i | cut -f 2 ; done | sort | uniq | wc -l`
  fastpNum=`ls ./*fastp*json | wc -l`
  mqcSampleNum=\$((samplesNum + fastpNum ))
  
  multiqc . --cl_config "max_table_rows: \$(( mqcSampleNum + 1 ))" -x ignoreFolder/ -x fastp_original/
  general_stats_parse.py --print-criteria 
  rm -rf multiqc_report.html multiqc_data

  if [ \$samplesNum -gt 50 ] ; then 
    cp genstats-QC_Status.txt QC_Status.txt
    beeswarm_config="max_table_rows: \${mqcSampleNum}"
  else
    beeswarm_config="max_table_rows: \$(( mqcSampleNum + 1 ))"
  fi

  multiqc . --cl_config "title: \\"Cohort MultiQC Report\\"" --cl_config "subtitle: \\"${cohort} QC\\"" --cl_config "intro_text: \\"Aggregate results from Tempo QC analysis\\"" --cl_config "\${beeswarm_config}" --cl_config "report_comment: \\"This report includes FASTQ and alignment for all samples in ${cohort}and Tumor/Normal pair statistics for all pairs in ${cohort}.\\"" -z -x ignoreFolder/ -x fastp_original/

  """
}

}
}

workflow.onComplete {
  file(params.fileTracking).text = ""
  file(outDir).eachFileRecurse{
    file(params.fileTracking).append(it + "\n")
  }
}

/*
================================================================================
=                          Reference Define Functions                          =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  if(params."${item}" == null){println "${item} is not found in reference map"; exit 1}
  if(params."${item}" == []){println "${item} is not found; glob pattern produces empty list"; exit 1}
  return file(params."${item}", checkIfExists: true)
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  result_array = [
    'dbsnp' : checkParamReturnFile("dbsnp"),
    'dbsnpIndex' : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict' : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile' : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex' : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex' : checkParamReturnFile("bwaIndex"), 
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels' : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
    'msiSensorList' : checkParamReturnFile("msiSensorList"),
    'svCallingExcludeRegions' : checkParamReturnFile("svCallingExcludeRegions"),
    'svCallingIncludeRegions' : checkParamReturnFile("svCallingIncludeRegions"),
    'svCallingIncludeRegionsIndex' : checkParamReturnFile("svCallingIncludeRegionsIndex"),
    // Target and Bait BED files
    'idtTargets' : checkParamReturnFile("idtTargets"),
    //'idtTargetsUnzipped' : checkParamReturnFile("idtTargetsUnzipped"),
    'idtTargetsIndex' : checkParamReturnFile("idtTargetsIndex"),
    'idtTargetsList' : checkParamReturnFile("idtTargetsList"),  
    'idtBaitsList' : checkParamReturnFile("idtBaitsList"), 
    'agilentTargets' : checkParamReturnFile("agilentTargets"),
    //'agilentTargetsUnzipped' : checkParamReturnFile("agilentTargetsUnzipped"),
    'agilentTargetsIndex' : checkParamReturnFile("agilentTargetsIndex"),
    'agilentTargetsList' : checkParamReturnFile("agilentTargetsList"),  
    'agilentBaitsList' : checkParamReturnFile("agilentBaitsList"), 
    'wgsTargets' : checkParamReturnFile("wgsTargets"),
    //'wgsTargetsUnzipped' : checkParamReturnFile("wgsTargetsUnzipped"),
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex")
  ]

  if (workflow.profile != "test") {
    result_array << ['vepCache' : checkParamReturnFile("vepCache")]
    // for SNP Pileup
    result_array << ['facetsVcf' : checkParamReturnFile("facetsVcf")]
    // intervals file for spread-and-gather processes
    result_array << ['intervals' : checkParamReturnFile("intervals")]
    // files for CombineChannel, needed by bcftools annotate
    result_array << ['repeatMasker' : checkParamReturnFile("repeatMasker")]
    result_array << ['repeatMaskerIndex' : checkParamReturnFile("repeatMaskerIndex")]
    result_array << ['mapabilityBlacklist' : checkParamReturnFile("mapabilityBlacklist")]
    result_array << ['mapabilityBlacklistIndex' : checkParamReturnFile("mapabilityBlacklistIndex")]
    // isoforms needed by vcf2maf
    result_array << ['isoforms' : checkParamReturnFile("isoforms")]
    // PON files
    result_array << ['exomePoN' : checkParamReturnFile("exomePoN")]
    result_array << ['exomePoNIndex' : checkParamReturnFile("exomePoNIndex")]
    result_array << ['wgsPoN' : checkParamReturnFile("wgsPoN")]
    result_array << ['wgsPoNIndex' : checkParamReturnFile("wgsPoNIndex")]
    // gnomAD resources
    result_array << ['gnomadWesVcf' : checkParamReturnFile("gnomadWesVcf")]
    result_array << ['gnomadWesVcfIndex' : checkParamReturnFile("gnomadWesVcfIndex")]
    result_array << ['gnomadWgsVcf' : checkParamReturnFile("gnomadWgsVcf")]
    result_array << ['gnomadWgsVcfIndex' : checkParamReturnFile("gnomadWgsVcfIndex")]
    // HLA FASTA and *dat for LOHHLA 
    result_array << ['hlaFasta' : checkParamReturnFile("hlaFasta")] 
    result_array << ['hlaDat' : checkParamReturnFile("hlaDat")] 
    // files for neoantigen & NetMHC
    result_array << ['neoantigenCDNA' : checkParamReturnFile("neoantigenCDNA")]
    result_array << ['neoantigenCDS' : checkParamReturnFile("neoantigenCDS")]
    // coding region BED files for calculating TMB
    result_array << ['idtCodingBed' : checkParamReturnFile("idtCodingBed")]
    result_array << ['agilentCodingBed' : checkParamReturnFile("agilentCodingBed")]    
    result_array << ['wgsCodingBed' : checkParamReturnFile("wgsCodingBed")]  
  }
  return result_array
}

def touchInputs() {
  new Timer().schedule({
  for ( i in epochMap.keySet() ){
    fileEpoch = file(i).lastModified()
    if ( fileEpoch > epochMap[i]) {
      epochMap[i] = fileEpoch
      "touch -ca ${i}".execute()
    }
  }
} as TimerTask, 1000, params.touchInputsInterval * 60 * 1000 ) // convert minutes to milliseconds
}

def watchMapping(tsvFile, assayType) {
  Channel.watchPath( tsvFile, 'create, modify' )
	 .splitCsv(sep: '\t', header: true)
	 .unique()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: false)
              def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: false)
              def numOfPairs = row.NUM_OF_PAIRS.toInteger()
              if(!TempoUtils.checkTarget(target, assayType)){}
              if(!TempoUtils.checkNumberOfItem(row, 5, tsvFile)){}

              [idSample, numOfPairs, target, fastqFile1, fastqFile2]
      }
      .map{ idSample, numOfPairs, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, numOfPairs), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchBamMapping(tsvFile, assayType){
  Channel.watchPath( tsvFile, 'create, modify' )
	 .splitCsv(sep: '\t', header: true)
	 .unique()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def bam = file(row.BAM, checkIfExists: false)
              def bai = file(row.BAI, checkIfExists: false)
              if(!TempoUtils.checkTarget(target, assayType)){}
              if(!TempoUtils.checkNumberOfItem(row, 4, tsvFile)){}

              [idSample, target, bam, bai]
      }
      .map{ idSample, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, 1), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchPairing(tsvFile){
  Channel.watchPath( tsvFile, 'create, modify' )
	 .splitCsv(sep: '\t', header: true)
	 .unique()
         .map { row ->
              def TUMOR_ID = row.TUMOR_ID
              def NORMAL_ID = row.NORMAL_ID
              if(!TempoUtils.checkNumberOfItem(row, 2, tsvFile)){}

              [TUMOR_ID, NORMAL_ID]
         }
	 .unique()
}

def watchAggregateWithPath(tsvFile) {
  Channel.watchPath(tsvFile, 'create, modify')
         .splitCsv(sep: '\t', header: true)
	 .unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              def path = row.PATH
              if(!TempoUtils.checkNumberOfItem(row, 5, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal, path]
         }
         .map { cohort, cohortSize, idTumor, idNormal, path
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal, path)
         }
         .transpose()
	 .unique()
}

def watchAggregate(tsvFile) {
  Channel.watchPath(file(runAggregate), 'create, modify')
         .splitCsv(sep: '\t', header: true)
	 .unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              if(!TempoUtils.checkNumberOfItem(row, 4, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal]
         }
         .map { cohort, cohortSize, idTumor, idNormal
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal)
         }
         .transpose()
	 .unique()
}
