#!/usr/bin/env nextflow

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

// Both mapping and pairing necessary for alignment of FASTQs
// Only bamPairing required when using already aligned BAM files
if (params.mapping && !params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --pairing and re-run the pipeline."
  exit 1
}

if (!params.mapping && params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --mapping and re-run the pipeline."
  exit 1
}

if ((params.mapping && params.bamPairing) || (params.pairing && params.bamPairing)) {
  println "ERROR: Cannot use both FASTQs and BAMs as inputs. Flags --bamPairing and --mapping/-pairing cannot be invoked together. Please provide either FASTQs or BAMs, and re-run the pipeline."
  exit 1
} 

// Validate mapping file
// Check for duplicate inputs and unique fileID
if (params.mapping) {
  mappingPath = params.mapping
  
  if (mappingPath && !TempoUtils.check_for_duplicated_rows(mappingPath)) {
    println "ERROR: Duplicated row found in mapping file. Please fix the error and re-run the pipeline."
    exit 1
  }

  // if (mappingPath && !TempoUtils.checkForUniqueSampleLanes(mappingPath)) {
  //   println "ERROR: The combination of sample ID and fileID values must be unique. Duplicate fileID for one sample cause errors. Please fix the error and re-run the pipeline."
  //   exit 1
  // }
}

// Validate pairing file
// Check for duplicate inputs
if (params.pairing) {
  pairingPath = params.pairing

  if (!TempoUtils.check_for_duplicated_rows(pairingPath)) {
    println "ERROR: Duplicated row found in pairing file. Please fix the error and re-run the pipeline."
    exit 1
  }
}

// Validate BAM file pairing file
// Check for duplicate inputs
if (params.bamPairing) {
  bamPairingPath = params.bamPairing

  if (bamPairingPath && !TempoUtils.check_for_duplicated_rows(bamPairingPath)) {
    println "ERROR: Duplicated row found in BAM mapping file. Please fix the error and re-run the pipeline."
    exit 1
  }
}

// User-set runtime parameters
publishAll = params.publishAll
outname = params.outname
runGermline = params.germline
runSomatic = params.somatic
runQC = params.QC

referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


// Skip these processes if starting from aligned BAM files
if (!params.bamPairing) {

  // Parse input FASTQ mapping and sample pairing
  mappingFile = file(mappingPath)
  pairingFile = file(pairingPath)
  TempoUtils.extractPairing(pairingFile).set{ inputPairing }
  TempoUtils.extractFastq(mappingFile).set{ inputFastqs }

  if (params.splitLanes) {
  inputFastqs
	.groupTuple(by: [0])
	.map{ idSample, fileID, files_pe1, files_pe2, targets
		-> tuple( groupKey(idSample, fileID.size()), fileID, files_pe1, files_pe2, targets)
	}
	.transpose()
	.into{ fastqsNeedSplit; fastqsNoNeedSplit }

  fastqsNeedSplit
	.filter{ item ->
		def idSample = item[0]
		def fileID = item[1]
		def file_pe1 = item[2]
		def file_pe2 = item[3]
		def targetFile = item[4]

		!(item[2].getName() =~ /_L(\d){3}_/)
	}
	.into{ inputFastqR1; inputFastqR2 }

  fastqsNoNeedSplit
	.filter{ item ->
		def idSample = item[0]
		def fileID = item[1]
		def file_pe1 = item[2]
		def file_pe2 = item[3]
		def targetFile = item[4]

		item[2].getName() =~ /_L(\d){3}_/
	}
	.map{ item ->
		def idSample = item[0]
		def fileID = item[1]
		def file_pe1 = item[2]
		def file_pe2 = item[3]
		def targetFile = item[4]

		return [ idSample, fileID, targetFile, file_pe1, file_pe2 ]
	}
	.set{ fastqsNoNeedSplit }

  process SplitLanesR1 {
    tag {idSample + "@" + fileID + "@R1"}   // The tag directive allows you to associate each process executions with a custom label

    input:
      set idSample, fileID, file(fastqFile1), file(fastqFile2), targetFile from inputFastqR1

    output:
      file("file-size.txt") into R1Size
      set idSample, fileID, file("*R1.splitLanes.fastq.gz"), targetFile into perLaneFastqsR1

    when: params.splitLanes

    script:
    inputSize = fastqFile1.size()
    if (workflow.profile == "juno") {
      if (inputSize > 10.GB) {
        task.time = { 72.h }
      }
      else if (inputSize < 5.GB) {
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
    }

    """
      echo -e "${idSample}@${fileID}@R1\t${inputSize}" > file-size.txt
      zcat $fastqFile1 | awk 'BEGIN {FS = ":"} {lane=\$4 ; print | "gzip > ${fileID}_L00"lane"_R1.splitLanes.fastq.gz" ; for (i = 1; i <= 3; i++) {getline ; print | "gzip > ${fileID}_L00"lane"_R1.splitLanes.fastq.gz"}}'
    """
  }
  process SplitLanesR2 {
    tag {idSample + "@" + fileID + "@R2"}   // The tag directive allows you to associate each process executions with a custom label

    input:
      set idSample, fileID, file(fastqFile1), file(fastqFile2), targetFile from inputFastqR2

    output:
      file("file-size.txt") into R2Size
      set idSample, fileID, file("*R2.splitLanes.fastq.gz"), targetFile into perLaneFastqsR2

    when: params.splitLanes

    script:
    inputSize = fastqFile2.size()
    if (workflow.profile == "juno") {
      if (inputSize > 10.GB) {
        task.time = { 72.h }
      }
      else if (inputSize < 5.GB) {
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
    }

    """
      echo -e "${idSample}@${fileID}@R2\t${inputSize}" > file-size.txt
      zcat $fastqFile2 | awk 'BEGIN {FS = ":"} {lane=\$4 ; print | "gzip > ${fileID}_L00"lane"_R2.splitLanes.fastq.gz" ; for (i = 1; i <= 3; i++) {getline ; print | "gzip > ${fileID}_L00"lane"_R2.splitLanes.fastq.gz"}}'
    """
  }

  fastqFiles = perLaneFastqsR1
	.combine(perLaneFastqsR2, by: [0,1,3])
	.concat(fastqsNoNeedSplit)
        .map{ item ->
            def idSample = item[0]
            def fileID = item[1]
            def file_pe1 = item[3]
            def file_pe2 = item[4]
            def targetFile = item[2]
	    def numOfLanes = file_pe1 instanceof Collection ? file_pe1.size() : 1

            return [ idSample, fileID, file_pe1, file_pe2, targetFile, numOfLanes ]
        }
	.groupTuple(by: [0])
	.map { idSample, fileID, file_pe1, file_pe2, targetFile, numOfLanes
		-> tuple(groupKey(idSample, numOfLanes.sum()), fileID, file_pe1, file_pe2, targetFile, numOfLanes)
	}
	.transpose()
        .transpose()
	.map{ item ->
	    def idSample = item[0]
	    def fileID = item[1] + "@" + TempoUtils.flowcellLaneFromFastq(item[2])[1]
	    def file_pe1 = item[2]
	    def file_pe1_size = item[2].size()
	    def file_pe2 = item[3]
	    def file_pe2_size = item[3].size()
	    def targetFile = item[4]
	    def rgID = TempoUtils.flowcellLaneFromFastq(item[2])[0] + ":" + TempoUtils.flowcellLaneFromFastq(item[2])[1]

	    return [ idSample, fileID, file_pe1, file_pe1_size, file_pe2, file_pe2_size, targetFile, rgID ]
	}
  }
  else{
     fastqFiles =  inputFastqs
        .map{ item ->
            def idSample = item[0]
            def fileID = item[1] + "@" + TempoUtils.flowcellLaneFromFastq(item[2])[1]
            def file_pe1 = item[2]
            def file_pe1_size = item[2].size()
            def file_pe2 = item[3]
            def file_pe2_size = item[3].size()
            def targetFile = item[4]
            def rgID = TempoUtils.flowcellLaneFromFastq(item[2])[0] + ":" + TempoUtils.flowcellLaneFromFastq(item[2])[1]

            return [ idSample, fileID, file_pe1, file_pe1_size, file_pe2, file_pe2_size, targetFile, rgID ]
        }
	.groupTuple(by:[0])
	.map{ idSample, fileID, files_pe1, files_pe1_size, files_pe2, files_pe2_size, targets, rgID
		-> tuple( groupKey(idSample, fileID.size()), fileID, files_pe1, files_pe1_size, files_pe2, files_pe2_size, targets, rgID)}
	.transpose()
  }


  // AlignReads - Map reads with BWA mem output SAM
  process AlignReads {
    tag {idSample + "@" + fileID}   // The tag directive allows you to associate each process executions with a custom label

    publishDir "${params.outDir}/qc/${idSample}/fastp", mode: params.publishDirMode, pattern: "*.{html,json}"

    input:
      set idSample, fileID, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, targetFile, rgID from fastqFiles
      set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

    output:
      file("*.html") into fastPHtml
      file("*.json") into fastPJson
      file("file-size.txt") into laneSize
      set idSample, fileID, file("${fileID}.sorted.bam"), targetFile into sortedBam

    script:
    // LSF resource allocation for juno
    // if running on juno, check the total size of the FASTQ pairs in order to allocate the runtime limit for the job, via LSF `bsub -W`
    // if total size of the FASTQ pairs is over 20 GB, use 72 hours
    // if total size of the FASTQ pairs is under 12 GB, use 3h. If there is a 140 error, try again with 6h. If 6h doesn't work, try 72h.
    inputSize = sizeFastqFile1 + sizeFastqFile2
    if (workflow.profile == "juno") {
      if (inputSize > 18.GB) {
        task.time = { 72.h }
      }
      else if (inputSize < 9.GB) {
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
    }

    // if it's the last time to try, use 72h as time limit no matter for what reason it failed before
    task.time = task.attempt < 3 ? task.time : { 72.h }
    
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

    readGroup = "@RG\\tID:${rgID}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
    """
    set -e
    set -o pipefail
    echo -e "${idSample}@${fileID}\t${inputSize}" > file-size.txt
    fastp --html ${fileID}.fastp.html --json ${fileID}.fastp.json --in1 ${fastqFile1} --in2 ${fastqFile2}
    bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} | samtools view -Sb - > ${fileID}.bam

    samtools sort -m ${mem}M -@ ${task.cpus} -o ${fileID}.sorted.bam ${fileID}.bam
    """
  }

  sortedBam.groupTuple()
    .map{ item ->
      def idSample = item[0]
      def fileID = item[1] //is a list
      def bam = item[2]
  
      def targetList = item[3].unique()

      if (targetList.size() > 1) {
        println "ERROR: Multiple targets found for ${idSample}; check inputs"
        exit 1
      }
    
      def target = targetList[0]

      [idSample, fileID, bam, target]
    }.set{ groupedBam }

  // MergeBams
  process MergeBams {
    tag {idSample}

    input:
      set idSample, fileID, file(bam), targetFile from groupedBam

    output:
      set idSample, file("${idSample}.merged.bam"), targetFile into mergedBam

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
    """
  }

  // GATK MarkDuplicates
  process MarkDuplicates {
    tag {idSample}

    input:
      set idSample, file(bam), targetFile from mergedBam

    output:
      set idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), targetFile into mdBams, mdBams4BQSR

    script:
    if (workflow.profile == "juno") {
      if(bam.size() > 120.GB) {
        task.time = { 72.h }
      }
      else if (bam.size() < 100.GB) {
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
    }
    // if it's the last time to try, use 72h as time limit no matter for what reason it failed before
    task.time = task.attempt < 3 ? task.time : { 72.h }

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
      set idSample, file(bam), file(bai), targetFile from mdBams
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
          task.time = { 72.h }
        }
        else if (bam.size() < 240.GB) {
          task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
        }
        else {
          task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
        }
      }
    }
    else {
      sparkConf = " BaseRecalibrator"
      task.cpus = 4
      task.memory = { 6.GB }
      task.time = { 72.h }
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

    publishDir "${params.outDir}/bams", mode: params.publishDirMode, pattern: "*.bam*"

    input:
      set idSample, file(bam), file(bai), target, file(recalibrationReport) from inputsBQSR
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict 
      ])

    output:
      set idSample, target, file("${idSample}.bam"), file("${idSample}.bam.bai") into bamsBQSR4Alfred, bamsBQSR4CollectHsMetrics, bamsBQSR4Tumor, bamsBQSR4Normal, bamsBQSR4QcPileup
      file("file-size.txt") into bamSize

    script:

    if (task.attempt < 3 ) {
      sparkConf = " ApplyBQSRSpark --conf 'spark.executor.cores = " + task.cpus + "'"
      if (workflow.profile == "juno") {
        if (bam.size() > 200.GB){
          task.time = { 72.h }
        }
        else if (bam.size() < 100.GB) {
          task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
        }
        else {
          task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
        }
      }
    }
    else {
      sparkConf = " ApplyBQSR"
      task.cpus = 4
      task.memory = { 6.GB }
      task.time = { 72.h }
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

  inputPairing.into{pairing4T; pairing4N; pairingTN}
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
			  .into{resultTsv; bamFiles}


  File file = new File(outname)
  file.newWriter().withWriter { w ->
      w << "TUMOR_ID\tNORMAL_ID\tTARGET\tTUMOR_BAM\tNORMAL_BAM\n"
  }

  if (workflow.profile == 'awsbatch') {
      resultTsv.subscribe { Object obj ->
        file.withWriterAppend { out ->
          out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\ts3:/${obj[3]}\ts3:/${obj[5]}"
        }
      }
    }
  else {
    resultTsv.subscribe { Object obj ->
      file.withWriterAppend { out ->
        out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\t${obj[5]}"
      }
    }
  }
  
}

/*
================================================================================
=                                SOMATIC PIPELINE                              =
================================================================================
*/

// parse --tools parameter for downstream 'when' conditionals, e.g. when: `` 'delly ' in tools
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
if ("conpair" in tools || params.conpairAll) {
  tools.add("pileup")
}

// If starting with BAM files, parse BAM pairing input
if (params.bamPairing) {
  bamPairingfile = file(bamPairingPath)
  TempoUtils.extractBAM(bamPairingfile).set { inputBams }
  inputBams.into{bamFiles; bamsTumor; bamsNormal; pairingTN}

  bamsTumor.map { item ->
		def idTumor = item[0]
		def idNormal = item[1]
		def tumorBam = item[3]
		def tumorBai = item[4]
		def target = item[2]
		return [ idTumor, idNormal, target, tumorBam, tumorBai ]
	}
	.unique()
	.into{bamsTumor4Combine; bamsTumor4VcfCombine}
  bamsNormal.map { item ->
		def idTumor = item[0]
		def idNormal = item[1]
		def normalBam = item[5]
		def normalBai = item[6]
		def target = item[2]
		 return [ idTumor, idNormal, target, normalBam, normalBai ]
	 }
	.unique()
	.into{ bamsNormal4Combine; bamsNormalOnly }
  bamsNormalOnly.map { item ->
			def idNormal = item[1]
			def target = item[2]
			def normalBam = item[3]
			def normalBai = item[4]
			return [ idNormal, target, normalBam, normalBai ] }
	.unique()
	.into{ bams4Haplotypecaller; bamsNormal4Polysolver; bamsForStrelkaGermline; bamsForMantaGermline; bamsForDellyGermline }


  bamsTumor4Combine.mix(bamsNormal4Combine)
		.map { item ->
			def target = item[2]
			def sampleBam = item[3]
			def sampleBai = item[4]
			def idSample = sampleBam.getSimpleName()
			return [ idSample, target, sampleBam, sampleBai ]
		}
		.into{bamsBQSR4Alfred; bamsBQSR4CollectHsMetrics; bamsBQSR4QcPileup}

  pairingTN.map{ item ->
		def idTumor = item[0]
		def idNormal = item[1]
		return [ idTumor, idNormal ]
	  }
	  .set { pairingTN }
}

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

bamFiles.into{bamsTN4Intervals; bamsForDelly; bamsForManta; bamsForMsiSensor; bamFiles4DoFacets; bamsForLOHHLA; }

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

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/mutect2", mode: params.publishDirMode

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

  input:
    each svType from svTypes
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}__${idNormal}_${svType}.filter.bcf") into dellyFilter4Combine

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
  """
}

// --- Run Manta
process SomaticRunManta {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/manta", mode: params.publishDirMode, pattern: "*.manta.vcf.{gz,gz.tbi}"

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
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

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

  publishDir "${params.outDir}/somatic/${outputPrefix}/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  publishDir "${params.outDir}/somatic/${outputPrefix}/merged_sv", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    set idTumor, idNormal, target, file(dellyBcfs), file(mantaFile) from dellyMantaCombineChannel

  output:
    set file("${outputPrefix}.delly.manta.vcf.gz"), file("${outputPrefix}.delly.manta.vcf.gz.tbi") into dellyMantaCombinedOutput, dellyMantaCombined4Aggregate
    set file("${outputPrefix}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz"), file("${outputPrefix}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz.tbi") into dellyOutput

  when: tools.containsAll(["manta", "delly"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """ 
  for f in *.bcf
  do 
    bcftools view --output-type z \$f > \${f%%.*}.delly.vcf.gz
  done

  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done
  
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

process SomaticRunStrelka2 {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/strelka2", mode: params.publishDirMode, pattern: "*.vcf.{gz,gz.tbi}"

  input:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])

  output:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file('*strelka2.vcf.gz'), file('*strelka2.vcf.gz.tbi') into strelkaAndBams4Combine
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


mutect2CombinedVcf4Combine.combine(strelkaAndBams4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }

// Combined Somatic VCFs

process SomaticCombineChannel {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutectCombinedVcf), file(mutectCombinedVcfIndex), file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(strelkaVcf), file(strelkaVcfIndex) from mutectStrelkaChannel
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

  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/maf_unfiltered/", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

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


mafFile.into{mafFileForMafAnno; mafFileForMutSig}

// --- Run Mutational Signatures, github.com/mskcc/mutation-signatures
process RunMutationSignatures {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.mutsig.txt") into mutSig4Aggregate

  when: tools.containsAll(["mutect2", "manta", "strelka2", "mutsig"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  python /mutation-signatures/main.py \
    /mutation-signatures/Stratton_signatures30.txt \
    ${outputPrefix}.somatic.maf \
    ${outputPrefix}.mutsig.txt
  """
}



// --- Run FACETS 
process DoFacets {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${tag}/facets", mode: params.publishDirMode, pattern: "*.snp_pileup.dat.gz"
  publishDir "${params.outDir}/somatic/${tag}/facets", mode: params.publishDirMode, pattern: "${tag}_OUT.txt"
  publishDir "${params.outDir}/somatic/${tag}/facets", mode: params.publishDirMode, pattern: "${outputDir}/*.{Rdata,png,seg,txt,out}"

  input:
    set idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal) from bamFiles4DoFacets
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set target, idTumor, idNormal, file("${outfile}") into snpPileupOutput
    set idTumor, idNormal, target, file("${outputDir}/*purity.out"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.seg"), file("${outputDir}/*hisens.out"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*purity.CNCF.png"), val("${outputFacetsSubdirectory}/${outputDir}") into FacetsOutput
    set idTumor, idNormal, target, file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*hisens.cncf.txt"), val("${outputDir}") into facetsForMafAnno, facetsForMafAnnoGermline
    set idTumor, idNormal, target, file("${outputDir}/*purity.out") into facetsPurity4LOHHLA, facetsPurity4MetaDataParser
    set file("${outputDir}/*purity.seg"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.CNCF.png"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.out") into FacetsPurity
    set file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.out") into FacetsHisens
    file("${tag}_OUT.txt") into FacetsPurityHisens4Aggregate

  when: "facets" in tools && runSomatic

  script:
  outfile = idTumor + "__" + idNormal + ".snp_pileup.dat.gz"
  tag = outputFacetsSubdirectory = "${idTumor}__${idNormal}"
  outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"
  """
  snp-pileup \
    --count-orphans \
    --pseudo-snps=50 \
    --gzip \
    ${facetsVcf} \
    ${outfile} \
    ${bamNormal} ${bamTumor}

  mkdir ${outputDir}

  Rscript --no-init-file /usr/bin/facets-suite/doFacets.R \
    --cval ${params.facets.cval} \
    --snp_nbhd ${params.facets.snp_nbhd} \
    --ndepth ${params.facets.ndepth} \
    --min_nhet ${params.facets.min_nhet} \
    --purity_cval ${params.facets.purity_cval} \
    --purity_snp_nbhd ${params.facets.purity_snp_nbhd} \
    --purity_ndepth ${params.facets.purity_ndepth} \
    --purity_min_nhet ${params.facets.purity_min_nhet} \
    --genome ${params.facets.genome} \
    --counts_file ${outfile} \
    --TAG ${tag} \
    --directory ${outputDir} \
    --R_lib /usr/lib/R/library \
    --seed ${params.facets.seed} \
    --tumor_id ${idTumor}

  python3 /usr/bin/facets-suite/summarize_project.py \
    -p ${tag} \
    -c ${outputDir}/*cncf.txt \
    -o ${outputDir}/*out \
    -s ${outputDir}/*seg  
  """
}


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
  tag {idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/lohhla", mode: params.publishDirMode

  input:
    set idNormal, target, idTumor, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal), file(purityOut), placeHolder, file(winnersHla) from mergedChannelLOHHLA
    set file(hlaFasta), file(hlaDat) from Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat])

  output:
    set file("*HLAlossPrediction_CI.txt"), file("*DNA.IntegerCPN_CI.txt"), file("*.pdf"), file("*.RData") optional true into lohhlaOutput

  when: tools.containsAll(["lohhla", "polysolver", "facets"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  cat ${winnersHla} | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
  
  PURITY=\$(grep Purity *_purity.out | grep -oP "[0-9\\.]+|NA+")
  PLOIDY=\$(grep Ploidy *_purity.out | grep -oP "[0-9\\.]+|NA+")
  cat <(echo -e "tumorPurity\ttumorPloidy") <(echo -e "\$PURITY\t\$PLOIDY") > tumor_purity_ploidy.txt

  Rscript --no-init-file /lohhla/LOHHLAscript.R \
    --patientId ${outputPrefix} \
    --normalBAMfile ${bamNormal} \
    --tumorBAMfile ${bamTumor} \
    --HLAfastaLoc ${hlaFasta} \
    --HLAexonLoc ${hlaDat} \
    --CopyNumLoc tumor_purity_ploidy.txt \
    --hlaPath massaged.winners.hla.txt \
    --gatkDir /picard-tools \
    --novoDir /opt/conda/bin

  if find Figures -mindepth 1 | read
  then
    mv Figures/* .
  fi
  """
}


facetsForMafAnno.combine(mafFileForMafAnno, by: [0,1,2]).set{ facetsMafFileSomatic }


// --- Do FACETS MAF annotation and post processing
process SomaticFacetsAnnotation {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/facets/${facetsPath}", mode: params.publishDirMode, pattern: "*{armlevel,genelevel}.unfiltered.txt"
  
  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileSomatic

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.facets.zygosity.maf") into FacetsAnnotationMafFile
    set idTumor, idNormal, target, file("${outputPrefix}.facets.zygosity.maf"), file("${outputPrefix}.armlevel.unfiltered.txt") into mafAndArmLevel4MetaDataParser
    set file("${outputPrefix}.armlevel.unfiltered.txt"), file("${outputPrefix}.genelevel.unfiltered.txt") into FacetsArmGene4Aggregate, FacetsArmGeneOutput
    file("file-size.txt") into mafSize

  when: tools.containsAll(["facets", "mutect2", "manta", "strelka2"]) && runSomatic

  script:
  mapFile = "${idTumor}__${idNormal}.map"
  outputPrefix = "${idTumor}__${idNormal}"
  """
  echo "Tumor_Sample_Barcode\tRdata_filename" > ${mapFile}
  echo "${idTumor}\t${purity_rdata.fileName}" >> ${mapFile}

  Rscript --no-init-file /usr/bin/facets-suite/mafAnno.R \
    --facets_files ${mapFile} \
    --maf ${maf} \
    --out_maf ${outputPrefix}.facets.maf
  
  Rscript --no-init-file /usr/bin/facets-suite/geneLevel.R \
    --filenames ${hisens_cncf} \
    --targetFile exome \
    --outfile ${outputPrefix}.genelevel.unfiltered.txt

  sed -i -e s@${idTumor}@${outputPrefix}@g ${outputPrefix}.genelevel.unfiltered.txt

  Rscript --no-init-file /usr/bin/facets-suite/armLevel.R \
    --filenames ${purity_cncf} \
    --outfile ${outputPrefix}.armlevel.unfiltered.txt

  sed -i -e s@${idTumor}@${outputPrefix}@g ${outputPrefix}.armlevel.unfiltered.txt

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-somatic.R ${outputPrefix}.facets.maf ${outputPrefix}.facets.zygosity.maf

  echo -e "${outputPrefix}\t`wc -l ${outputPrefix}.facets.zygosity.maf | cut -d ' ' -f1`" > file-size.txt

  """
}


hlaOutput.combine(FacetsAnnotationMafFile, by: [1,2]).set{ input4Neoantigen }

// --- Run neoantigen prediction pipeline
process RunNeoantigen {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/maf_final/", mode: params.publishDirMode, pattern: "*.maf"
  publishDir "${params.outDir}/somatic/${outputPrefix}/neoantigen/", mode: params.publishDirMode, pattern: "*.txt"

  input:
    set idNormal, target, placeHolder, file(polysolverFile), idTumor, file(mafFile) from input4Neoantigen
    set file(neoantigenCDNA), file(neoantigenCDS) from Channel.value([
      referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS
    ])

  output:
    file("${idTumor}__${idNormal}.all_neoantigen_predictions.txt") into NetMhcStats4Aggregate
    file("${idTumor}__${idNormal}.all_neoantigen_predictions.txt") into NetMhcStatsOutput
    file("*.maf") into NeoantigenMaf4Aggregate
    file("*.maf") into NeoantigenMafOutput

  when: tools.containsAll(["neoantigen", "mutect2", "manta", "strelka2"]) && runSomatic

  script:

  if (workflow.profile == "juno") {
    if(mafFile.size() > 10.MB){
      task.time = { 72.h }
    }
    else if (mafFile.size() < 5.MB){
      task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
    }
    else {
      task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
    }
  }

  outputPrefix = "${idTumor}__${idNormal}"
  tmpDir = "neoantigen-tmp"
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
    --output_dir ./

  awk 'NR==1 {printf("%s\\t%s\\n", "sample", \$0)} NR>1 {printf("%s\\t%s\\n", "${outputPrefix}", \$0) }' *.all_neoantigen_predictions.txt > ${outputPrefix}.all_neoantigen_predictions.txt
  mv ${outputPrefix}.neoantigens.maf ${outputPrefix}.somatic.final.maf
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
    set idTumor, idNormal, target, file("${outputPrefix}.msisensor.tsv") into msi4MetaData

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


facetsPurity4MetaDataParser.combine(mafAndArmLevel4MetaDataParser, by: [0,1,2])
			   .combine(msi4MetaData, by: [0,1,2])
			   .combine(mutSig4Aggregate, by: [0,1,2])
			   .combine(hlaOutputForMetaDataParser, by: [1,2])
			   .unique()
			   .set{ mergedChannelMetaDataParser }

// --- Generate sample-level metadata
process MetaDataParser {
  tag {idTumor + "__" + idNormal}
 
  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/meta_data/", mode: params.publishDirMode, pattern: "*.sample_data.txt"

  input:
    set idNormal, target, idTumor, file(purityOut), file(mafFile), file(armLevel), file(msifile), file(mutSig), placeHolder, file(polysolverFile) from mergedChannelMetaDataParser
    set file(idtCodingBed), file(agilentCodingBed), file(wgsCodingBed) from Channel.value([
      referenceMap.idtCodingBed, referenceMap.agilentCodingBed, referenceMap.wgsCodingBed
    ]) 

  output:
    file("*.sample_data.txt") into MetaDataOutput
    file("*.sample_data.txt") into MetaData4Aggregate

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
    --facetsArmLevel ${armLevel} \
    --MSIsensor_output ${msifile} \
    --mutational_signatures_output ${mutSig} \
    --polysolver_output ${polysolverFile} \
    --MAF_input ${mafFile} \
    --coding_baits_BED ${codingRegionsBed}
  
  mv ${idTumor}__${idNormal}_metadata.txt ${idTumor}__${idNormal}.sample_data.txt
  """
}


/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/

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

  publishDir "${params.outDir}/germline/${idNormal}/haplotypecaller", mode: params.publishDirMode

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

  publishDir "${params.outDir}/germline/${idNormal}/strelka2", mode: params.publishDirMode

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

  publishDir "${params.outDir}/germline/${idNormal}/maf_unfiltered", mode: params.publishDirMode, pattern: "*.unfiltered.maf"

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

  publishDir "${params.outDir}/germline/${idNormal}/maf_final/", mode: params.publishDirMode, pattern: "*.germline.final.maf"

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileGermline

  output:
    file("${outputPrefix}.germline.final.maf") into mafFileOutputGermline
    file("${outputPrefix}.germline.final.maf") into mafFile4AggregateGermline

  when: tools.containsAll(["facets", "haplotypecaller", "strelka2"]) && runGermline

  script:
  mapFile = "${idTumor}_${idNormal}.map"
  outputPrefix = "${idTumor}__${idNormal}.germline"
  """
  echo "Tumor_Sample_Barcode\tRdata_filename" > ${mapFile}
  echo "${idTumor}\t${purity_rdata.fileName}" >> ${mapFile}

  Rscript --no-init-file /usr/bin/facets-suite/mafAnno.R \
    --facets_files ${mapFile} \
    --maf ${maf} \
    --out_maf ${outputPrefix}.facets.maf

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-germline.R ${outputPrefix}.facets.maf ${outputPrefix}.germline.final.maf
  """
}

// --- Call germline SVs with Delly
Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypesGermline }

process GermlineDellyCall {
  tag {idNormal + '@' + svType}

  input:
    each svType from svTypesGermline
    set idNormal, target, file(bamNormal), file(baiNormal) from bamsForDellyGermline
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idNormal, target, file("${idNormal}_${svType}.filter.bcf") into dellyFilter4CombineGermline

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
  """
}

// --- Run Manta, germline
process GermlineRunManta {
  tag {idNormal}

  publishDir "${params.outDir}/germline/${idNormal}/manta", mode: params.publishDirMode

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

  publishDir "${params.outDir}/germline/${idNormal}/delly", mode: params.publishDirMode, pattern: "*delly.vcf.{gz,gz.tbi}"
  publishDir "${params.outDir}/germline/${idNormal}/combined_sv/", mode: params.publishDirMode, pattern: "*.delly.manta.vcf.{gz,gz.tbi}"

  input:
    set idNormal, target, file(dellyBcf), file(mantaVcf), file(mantaVcfIndex) from dellyMantaChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set file("${idNormal}.delly.manta.vcf.gz"), file("${idNormal}.delly.manta.vcf.gz.tbi") into dellyMantaCombinedOutputGermline, dellyMantaCombined4AggregateGermline
    set file("*delly.vcf.gz"), file("*delly.vcf.gz.tbi") into dellyOutputGermline

  when: tools.containsAll(["manta", "delly"]) && runGermline

  script:
  """ 
  for f in *.bcf
  do 
    bcftools view --output-type z \$f > \${f%%.*}.delly.vcf.gz
  done

  for f in *.delly.vcf.gz
  do
    tabix --preset vcf \$f
  done

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



/*
================================================================================
=                              Quality Control                                 =
================================================================================
*/

// GATK CollectHsMetrics, WES only
process QcCollectHsMetrics {
  tag {idSample}

  publishDir "${params.outDir}/qc/${idSample}/collecthsmetrics", mode: params.publishDirMode

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
    file("${idSample}.hs_metrics.txt") into collectHsMetricsOutput
    file("${idSample}.hs_metrics.txt") into collectHsMetrics4Aggregate

  when: params.assayType == "exome" && !params.test && runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { 72.h }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
    }
    else {
      task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
    }
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

// Alfred, BAM QC
Channel.from(true, false).set{ ignore_read_groups }
process QcAlfred {
  tag {idSample + "@" + "ignore_rg_" + ignore_rg }

  publishDir "${params.outDir}/qc/${idSample}/alfred", mode: params.publishDirMode

  input:
    each ignore_rg from ignore_read_groups
    set idSample, target, file(bam), file(bai) from bamsBQSR4Alfred
    file(genomeFile) from Channel.value([referenceMap.genomeFile])
    set file(idtTargets), file(agilentTargets), file(idtTargetsIndex), file(agilentTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex
    ])

  output:
    file("${idSample}.alfred*tsv.gz") into bamsQcStats4Aggregate
    set file("${idSample}.alfred*tsv.gz"), file("${idSample}.alfred*tsv.gz.pdf") into alfredOutput

  when: runQC

  script:
  if (workflow.profile == "juno") {
    if (bam.size() > 200.GB) {
      task.time = { 72.h }
    }
    else if (bam.size() < 100.GB) {
      task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
    }
    else {
      task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
    }
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


process QcPileup {
  tag {idSample}

  publishDir "${params.outDir}/qc/${idSample}/pileup/", mode: params.publishDirMode

  input:
    set idSample, target, file(bam), file(bai) from bamsBQSR4QcPileup
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idSample, file("${idSample}.pileup") into pileupOutput, tumorPileups, normalPileups

  when: !params.test && "pileup" in tools && runQC

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


pileupT.combine(pileupN, by: [0, 1]).set{ pileupConpair }

process QcConpair {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outPrefix}/conpair/", mode: params.publishDirMode

  input:
    set idTumor, idNormal, file(pileupTumor), file(pileupNormal) from pileupConpair
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    file("${outPrefix}.concordance.txt") into conpairConcordance
    file("${outPrefix}.contamination.txt") into conpairContamination
    set file("${outPrefix}.concordance.txt"), file("${outPrefix}.contamination.txt") into conpairOutput

  when: !params.test && "conpair" in tools && runQC

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

pileupT4Combine.combine(pileupN4Combine).set{ pileupConpairAll }

process QcConpairAll {
  tag {idTumor + "@" + idNormal}

  input:
    set idTumor, idNormal_noUse, file(pileupTumor), idTumor_noUse, idNormal, file(pileupNormal) from pileupConpairAll
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    file("${outPrefix}.concordance.txt") into conpairAllConcordance
    file("${outPrefix}.contamination.txt") into conpairAllContamination
    set file("${outPrefix}.concordance.txt"), file("${outPrefix}.contamination.txt") into conpairAllOutput

  when: !params.test && params.conpairAll && runQC

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


/*
================================================================================
=                              Cohort Aggregation                              =
================================================================================
*/

process SomaticAggregateMaf {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(mafFile) from NeoantigenMaf4Aggregate.collect()

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
  cat mut/*.maf | grep ^Hugo_Symbol | head -n 1 > mut_somatic.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo_Symbol" | sort -k5,5V -k6,6n >> mut_somatic.maf
  """
}

process SomaticAggregateNetMHC {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(netmhcCombinedFile) from NetMhcStats4Aggregate.collect()

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

process SomaticAggregateFacets {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(purityFiles) from FacetsPurity.collect()
    file(hisensFiles) from FacetsHisens.collect()
    file(purityHisens) from FacetsPurityHisens4Aggregate.collect()
    file(annotationFiles) from FacetsArmGene4Aggregate.collect()

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
  mv *{genelevel,armlevel}.unfiltered.txt facets_tmp/
  cat facets_tmp/*genelevel.unfiltered.txt | head -n 1 > cna_genelevel.txt
  awk -v FS='\t' '{ if (\$16 != "DIPLOID" && (\$17 == "PASS" || (\$17 == "FAIL" && \$18 == "rescue")))  print \$0 }' facets_tmp/*genelevel.unfiltered.txt >> cna_genelevel.txt
  cat facets_tmp/*armlevel.unfiltered.txt | head -n 1 > cna_armlevel.txt
  cat facets_tmp/*armlevel.unfiltered.txt | grep -v "DIPLOID" | grep -v "Tumor_Sample_Barcode" >> cna_armlevel.txt || [[ \$? == 1 ]]
  """
}

process SomaticAggregateSv {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(dellyMantaVcf) from dellyMantaCombined4Aggregate.collect()

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

process SomaticAggregateMetadata {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(metaDataFile) from MetaData4Aggregate.collect()

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


// --- Aggregate per-sample germline data, MAF
process GermlineAggregateMaf {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(mafFile) from mafFile4AggregateGermline.collect()

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
  cat mut/*.maf | grep ^Hugo | head -n1 > mut_germline.maf
  cat mut/*.maf | grep -Ev "^#|^Hugo" | sort -k5,5V -k6,6n >> mut_germline.maf

  """
}


// --- Aggregate per-sample germline data, SVs
process GermlineAggregateSv {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(dellyMantaVcf) from dellyMantaCombined4AggregateGermline.collect()

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


process QcBamAggregate {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(metricsFile) from collectHsMetrics4Aggregate.collect()
    file(bamsQcStatsFile) from bamsQcStats4Aggregate.collect()

  output:
    file('alignment_qc.txt') into alignmentQcAggregatedOutput

  when: !params.test && runQC

  script:
  if (params.assayType == "exome") {
    options = "wes"
  }
  else {
    options = 'wgs'
  }
  """
  Rscript --no-init-file /usr/bin/create-aggregate-qc-file.R ${options}
  """
}

// -- Run based on QcConpairAll channels or the single QcConpair channels
(conpairAggregateConcordance, conpairAggregateContamination) = (!params.conpairAll
								? [conpairConcordance, conpairContamination]
								: [conpairAllConcordance, conpairAllContamination]
								)

process QcConpairAggregate {

  publishDir "${params.outDir}/cohort_level", mode: params.publishDirMode

  input:
    file(concordance) from conpairAggregateConcordance.collect()
    file(contamination) from conpairAggregateContamination.collect()

  output:
    set file('concordance_qc.txt'), file('contamination_qc.txt') into conpairAggregatedOutput

  when: !params.test && runQC

  script:
  """
  echo -e "Pair\tConcordance" > concordance_qc.txt
  grep -v "concordance" *.concordance.txt | sed 's/.concordance.txt:/\t/' | cut -f1,3 | sort -k1,1 >> concordance_qc.txt
  echo -e "Pair\tSample_Type\tSample_ID\tContamination" > contamination_qc.txt
  grep -v "Contamination" *.contamination.txt | sed 's/.contamination.txt:/\t/' | sort -k1,1 >> contamination_qc.txt
  """
}



/*
================================================================================
=                          Reference Define Functions                          =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
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
    'idtTargetsIndex' : checkParamReturnFile("idtTargetsIndex"),
    'idtTargetsList' : checkParamReturnFile("idtTargetsList"),  
    'idtBaitsList' : checkParamReturnFile("idtBaitsList"), 
    'agilentTargets' : checkParamReturnFile("agilentTargets"),
    'agilentTargetsIndex' : checkParamReturnFile("agilentTargetsIndex"),
    'agilentTargetsList' : checkParamReturnFile("agilentTargetsList"),  
    'agilentBaitsList' : checkParamReturnFile("agilentBaitsList"), 
    'wgsTargets' : checkParamReturnFile("wgsTargets"),
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex")
  ]

  if (!params.test) {
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
