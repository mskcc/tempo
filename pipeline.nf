#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
Processes overview:

Alignment and QC
----------------
 - AlignReads
    --- Map paired-end FASTQs with bwa mem
    --- Sort BAM with samtools sort
    --- FASTQ QC with FastP 
 - MergeBam --- Merge BAM for the same samples from different lanes, samtools merge
 - MarkDuplicates --- Mark Duplicates with GATK4 MarkDuplicates
 - CreateRecalibrationTable --- Create Recalibration Table with GATK4 BaseRecalibrator
 - RecalibrateBam --- Recalibrate Bam with GATK4 ApplyBQSR
 - Alfred - BAM QC metrics
 - CollectHsMetrics --- *For WES only* Calculate hybrid-selection metrics, GATK4 CollectHsMetrics
 - AggregateBamQC --- aggregates information from Alfred and CollectHsMetrics across all samples

Somatic Analysis
----------------
 - CreateScatteredIntervals --- GATK4 SplitIntervals
 - RunMutect2 --- somatic SNV calling, MuTect2
 - SomaticRunStrelka2 --- somatic SNV calling, Strelka2, using Manta for small indel calling by default
 - SomaticCombineMutect2Vcf --- combine Mutect2 calls, bcftools
 - SomaticRunManta --- somatic SV calling, Manta
 - SomaticDellyCall --- somatic SV calling, Delly
 - SomaticMergeDellyAndManta --- combine Manta and Delly VCFs
 - SomaticCombineChannel --- combine and filter VCFs, bcftools
 - SomaticAnnotateMaf --- annotate MAF, vcf2maf
 - RunMsiSensor --- MSIsensor
 - DoFacets --- facets-suite: mafAnno.R, geneLevel.R, armLevel.R
 - RunPolysolver --- Polysolver
 - RunLOHHLA --- LOH in HLA
 - RunConpair --- Tumor-Normal quality/contamination
 - RunMutationSignatures --- mutational signatures
 - SomaticFacetsAnnotation --- annotate FACETS
 - RunNeoantigen --- NetMHCpan 4.0
 - MetaDataParser --- python script to parse metadata into single *tsv
 - SomaticAggregateMaf --- collect outputs, MAF
 - SomaticAggregateNetMHC --- collect outputs, neoantigen prediction
 - SomaticAggregateFacets --- collect outputs, FACETS
 - SomaticAggregateSv --- collect outputs, SVs
 - SomaticAggregateMetaData --- collect outputs, sample data

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
 - GermlineAggregateMaf --- collect outputs, MAF
 - GermlineAggregateSv --- collect outputs, SVs

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
// Only bam_pairing required when using already aligned BAM files
if (params.mapping && !params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --pairing and re-run the pipeline."
  exit 1
}

if (!params.mapping && params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --mapping and re-run the pipeline."
  exit 1
}

if ((params.mapping && params.bam_pairing) || (params.pairing && params.bam_pairing)) {
  println "ERROR: Cannot use both FASTQs and BAMs as inputs. Flags --bam_pairing and --mapping/-pairing cannot be invoked together. Please provide either FASTQs or BAMs, and re-run the pipeline."
  exit 1
} 

// Validate mapping file
// Check for duplicate inputs, mixed assay types and unique lane names
if (params.mapping) {
  mappingPath = params.mapping
  
  if (mappingPath && !TempoUtils.check_for_duplicated_rows(mappingPath)) {
    println "ERROR: Duplicated row found in mapping file. Please fix the error and re-run the pipeline."
    exit 1
  }

  if (mappingPath && !TempoUtils.check_for_mixed_assay(mappingPath)) {
    println "ERROR: Multiple assays found in mapping file. Users can either run exomes or genomes, but not both. Please fix the error and re-run the pipeline."
    exit 1
  }

  if (mappingPath && !TempoUtils.checkForUniqueSampleLanes(mappingPath)) {
    println "ERROR: The combination of sample ID and lane names values must be unique. Duplicate lane names for one sample cause errors. Please fix the error and re-run the pipeline."
    exit 1
  }
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
if (params.bam_pairing) {
  bamPairingPath = params.bam_pairing

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

referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


// Skip these processes if starting from aligned BAM files
if (!params.bam_pairing) {

  // Parse input FASTQ mapping and sample pairing
  fastqFiles = Channel.empty() 
  mappingFile = file(mappingPath)
  pairingFile = file(pairingPath)
  pairingTN = TempoUtils.extractPairing(pairingFile)
  fastqFiles = TempoUtils.extractFastq(mappingFile)

  fastqFiles =  fastqFiles.groupTuple(by:[0]).map{ key, lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets -> tuple( groupKey(key, lanes.size()), lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets)}.transpose()

  // AlignReads - Map reads with BWA mem output SAM
  process AlignReads {
    tag {idSample + "@" + lane}   // The tag directive allows you to associate each process executions with a custom label

    publishDir "${params.outDir}/qc/fastp/${idSample}", mode: params.publishDirMode, pattern: "*.html"
    if (publishAll) { 
      publishDir "${params.outDir}/qc/fastp/json", mode: params.publishDirMode, pattern: "*.json" 
    }

    input:
      set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assay, targetFile from fastqFiles
      set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

    output:
      file("*.html") into fastPHtml
      file("*.json") into fastPJson
      set idSample, lane, file("${lane}.sorted.bam"), assay, targetFile into sortedBam

    script:

// LSF resource allocation for juno
// if running on juno, check the total size of the FASTQ pairs in order to allocate the runtime limit for the job, via LSF `bsub -W`
// if total size of the FASTQ pairs is over 20 GB, use 72 hours
// if total size of the FASTQ pairs is under 12 GB, use 3h. If there is a 140 error, try again with 6h. If 6h doesn't work, try 72h.

    inputSize = sizeFastqFile1 + sizeFastqFile2
    if (workflow.profile == "juno") {
      if(inputSize > 20.GB){
        task.time = { 72.h }
      }
      else if (inputSize < 12.GB){
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
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

    readGroup = "@RG\\tID:${lane}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
    """
    set -e
    set -o pipefail
    echo -e "${lane}\t${inputSize}" > file-size.txt
    fastp --html ${lane}.fastp.html --in1 ${fastqFile1} --in2 ${fastqFile2}
    bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} | samtools view -Sb - > ${lane}.bam

    samtools sort -m ${mem}M -@ ${task.cpus} -o ${lane}.sorted.bam ${lane}.bam
    """
  }


  sortedBam.groupTuple().set{ groupedBam }

  groupedBam = groupedBam.map{ item -> 
    def idSample = item[0]
    def lane = item[1] //is a list
    def bam = item[2]
  
    def assayList = item[3].unique()
    def targetList = item[4].unique()

    if ((assayList.size() > 1) || (targetList.size() > 1)) {  
      println "ERROR: Multiple assays and/or targets found for ${idSample}; check inputs"
      exit 1
    }
    
    def assay = assayList[0]
    def target = targetList[0]

    [idSample, lane, bam, assay, target]
  }

  // MergeBams
  process MergeBams {
    tag {idSample}

    input:
      set idSample, lane, file(bam), assay, targetFile from groupedBam

    output:
      set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile into mergedBam

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
    """
  }

  // GATK MarkDuplicates
  process MarkDuplicates {
    tag {idSample}

    input:
      set idSample, lane, file(bam), assay, targetFile from mergedBam

    output:
      set file("${idSample}.md.bam"), file("${idSample}.md.bai"), idSample, lane, assay, targetFile into duplicateMarkedBams
      set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), assay, targetFile into markDuplicatesTSV
      file ("${idSample}.bam.metrics") into markDuplicatesReport

    script:
    if (workflow.profile == "juno" && params.assayType == "exome") {
      if(bam.size() > 200.GB) {
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
    maxMem = (memMultiplier * task.memory.toString().split(" ")[0].toInteger() - 3)
    maxMem = maxMem < 4 ? 5 : maxMem
    javaOptions = "--java-options '-Xms4000m -Xmx" + maxMem + "g'"
    """
    gatk MarkDuplicates \
      ${javaOptions} \
      --TMP_DIR ${TMPDIR} \
      --MAX_RECORDS_IN_RAM 50000 \
      --INPUT ${idSample}.merged.bam \
      --METRICS_FILE ${idSample}.bam.metrics \
      --ASSUME_SORT_ORDER coordinate \
      --CREATE_INDEX true \
      --OUTPUT ${idSample}.md.bam
    """
  }

  duplicateMarkedBams = duplicateMarkedBams.map {
      bam, bai, idSample, lane, assay, targetFile -> tag = bam.baseName.tokenize('.')[0]
      [idSample, bam, bai, assay, targetFile]
  }

  (mdBam, mdBamToJoin) = duplicateMarkedBams.into(2)

 // GATK BaseRecalibrator , CreateRecalibrationTable 
  process CreateRecalibrationTable {
    tag {idSample}

    input:
      set idSample, file(bam), file(bai), assay, targetFile from mdBam 
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
      set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), val("${idSample}.recal.table"), assay, targetFile into recalibrationTableTSV

    script:
    if (task.attempt < 3){
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
      task.cpus = 1
      task.memory = { 4.GB }
      task.time = { 72.h }
    }

    memMultiplier = params.mem_per_core ? task.cpus : 1
    javaOptions = "--java-options '-Xmx" + task.memory.toString().split(" ")[0].toInteger() * memMultiplier + "g'"
    knownSites = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
    """
    gatk \
      ${sparkConf} \
      ${javaOptions} \
      --tmp-dir ${TMPDIR} \
      --reference ${genomeFile} \
      --known-sites ${dbsnp} \
      ${knownSites} \
      --verbosity INFO \
      --input ${bam} \
      --output ${idSample}.recal.table
    """ 
  }

  recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0])

  // GATK ApplyBQSR, RecalibrateBAM
  process RecalibrateBam {
    tag {idSample}

    publishDir "${params.outDir}/bams", mode: params.publishDirMode

    input:
      set idSample, file(bam), file(bai), assay, targetFile, file(recalibrationReport) from recalibrationTable
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict 
      ])

    output:
      set idSample, file("${idSample}.bam"), file("${idSample}.bam.bai"), assay, targetFile into recalibratedBam, recalibratedBamForCollectHsMetrics, recalibratedBamForStats, recalibratedBamForOutput, recalibratedBamForOutput2
      file("${idSample}.bam") into currentBam
      file("${idSample}.bam.bai") into currentBai
      val(assay) into assays
      val(targetFile) into targets

    script:

    if (task.attempt < 3){
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
      task.cpus = 1
      task.memory = { 4.GB }
      task.time = { 72.h }
    }
    memMultiplier = params.mem_per_core ? task.cpus : 1
    javaOptions = "--java-options '-Xmx" + task.memory.toString().split(" ")[0].toInteger() * memMultiplier + "g'"
    """
    echo -e "${idSample}\t${bam.size()}" > file-size.txt
    gatk \
      ${sparkConf} \
      ${javaOptions} \
      --tmp-dir ${TMPDIR} \
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

  // set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
  recalibratedBamForOutput.combine(pairingTN)
                          .filter { item -> // only keep combinations where sample is same as tumor pair sample
                            def idSample = item[0]
                            def sampleBam = item[1]
                            def sampleBai = item[2]
                            def assay = item[3]
                            def target = item[4]
                            def idTumor = item[5]
                            def idNormal = item[6]
                            idSample == idTumor
                          }.map { item -> // re-order the elements
                            def idSample = item[0]
                            def sampleBam = item[1]
                            def sampleBai = item[2]
                            def assay = item[3]
                            def target = item[4]
                            def idTumor = item[5]
                            def idNormal = item[6]
                            def bamTumor = sampleBam
                            def baiTumor = sampleBai

                            return [ assay, target, idTumor, idNormal, bamTumor, baiTumor ]
                          }.combine(recalibratedBamForOutput2)
                          .filter { item ->
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def baiTumor = item[5]
                            def idSample = item[6]
                            def bamNormal = item[7]
                            def baiNormal = item[8]
                            idSample == idNormal
                          }.map { item -> // re-order the elements
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def baiTumor = item[5]
                            def idSample = item[6]
                            def bamNormal = item[7]
                            def baiNormal = item[8]

                            return [ assay, target, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
                          }
                          .set { result }

  result.into { resultTsv; bamFiles }

  File file = new File(outname)
  file.newWriter().withWriter { w ->
      w << "ASSAY\tTARGET\tTUMOR_ID\tNORMAL_ID\tTUMOR_BAM\tNORMAL_BAM\tTUMOR_BAI\tNORMAL_BAI\n"
  }

  if (workflow.profile == 'awsbatch') {
      resultTsv.subscribe { Object obj ->
        file.withWriterAppend { out ->
          out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\ts3:/${obj[4]}\ts3:/${obj[5]}\ts3:/${obj[6]}\ts3:/${obj[7]}"
        }
      }
    }
  else {
    resultTsv.subscribe { Object obj ->
      file.withWriterAppend { out ->
        out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\t${obj[4]}\t${obj[5]}\t${obj[6]}\t${obj[7]}"
      }
    }
  }
  
  // GATK CollectHsMetrics, WES only
  process CollectHsMetrics {
    tag {idSample}

    publishDir "${params.outDir}/qc/collecthsmetrics/${idSample}", mode: params.publishDirMode

    input:
      set idSample, file(bam), file(bai), assay, target from recalibratedBamForCollectHsMetrics
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
      ])
      set file(idtTargetsList), file(agilentTargetsList), file(idtBaitsList), file(agilentBaitsList) from Channel.value([
        referenceMap.idtTargetsList, referenceMap.agilentTargetsList, 
        referenceMap.idtBaitsList, referenceMap.agilentBaitsList      
      ])

    output:
      file("${idSample}.hs_metrics.txt") into collectHsMetrics

    when: 'wes' in assay && !params.test

    script:
    if (workflow.profile == "juno") {
      if(bam.size() > 200.GB){
        task.time = { 72.h }
      }
      else if (bam.size() < 100.GB){
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
      --TMP_DIR ${TMPDIR} \
      --INPUT ${bam} \
      --OUTPUT ${idSample}.hs_metrics.txt \
      --REFERENCE_SEQUENCE ${genomeFile} \
      --BAIT_INTERVALS ${baitIntervals} \
      --TARGET_INTERVALS ${targetIntervals} 
    """
  }

  // Alfred, BAM QC
  ignore_read_groups = Channel.from(true, false)
  process Alfred {
    tag {idSample + "@" + "ignore_rg_" + ignore_rg }

    publishDir "${params.outDir}/qc/alfred/${idSample}", mode: params.publishDirMode
  
    input:
      each ignore_rg from ignore_read_groups
      set idSample, file(bam), file(bai), assay, target from recalibratedBam
      file(genomeFile) from Channel.value([referenceMap.genomeFile])
      set file(idtTargets), file(agilentTargets), file(idtTargetsIndex), file(agilentTargetsIndex) from Channel.value([
        referenceMap.idtTargets, referenceMap.agilentTargets,
        referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex
      ])

    output:
      file("${idSample}.alfred*tsv.gz") into bamsQcStats
      file("${idSample}.alfred*tsv.gz.pdf") into bamsQcPdfs

    script:
    if (workflow.profile == "juno" && params.assayType == "exome") {
      if(bam.size() > 200.GB){
        task.time = { 72.h }
      }
      else if (bam.size() < 100.GB){
        task.time = task.exitStatus != 140 ? { 3.h } : { 6.h }
      }
      else {
        task.time = task.exitStatus != 140 ? { 6.h } : { 72.h }
      }
    }

    options = ""
    if (assay == "wes") {
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
  
  process AggregateBamQc {
    
    publishDir "${params.outDir}/qc", mode: params.publishDirMode

    input:
      file(metricsFile) from collectHsMetrics.collect()
      file(bamsQcStatsFile) from bamsQcStats.collect()

    output:
      file('alignment_qc.txt') into alignmentQc

    when: !params.test

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
if ("conpair" in tools || "conpairAll" in tools) {
  tools.add("pileup")
}

// If starting with BAM files, parse BAM pairing input
if (params.bam_pairing) {
  bamFiles = Channel.empty()
  bamPairingfile = file(bamPairingPath)
  bamFiles = TempoUtils.extractBAM(bamPairingfile)
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
    set file("agilent*.interval_list"), val("agilent") into agilentIList
    set file("idt*.interval_list"), val("idt") into idtIList
    set file("wgs*.interval_list"), val("wgs") into wgsIList

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


(bamsForIntervals, bamFiles) = bamFiles.into(2)


//Associating interval_list files with BAM files, putting them into one channel

(aBamList, iBamList, wBamList) = bamsForIntervals.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 1)
iMergedChannel = iBamList.combine(idtIList, by: 1)
wMergedChannel = wBamList.combine(wgsIList, by: 1)

// These will go into mutect2 and haplotypecaller

// From Nextflow Doc: However there are use cases in which each tuple has a different size depending grouping key. In this cases use the built-in function groupKey that allows you to create a special grouping key object to which it's possible to associate the group size for a given key.
// Reference: https://github.com/nextflow-io/nextflow/issues/796
// using groupKey() function to create a unique key for each TN pair, and the key also includes number of intervalBeds for each samples
// this change will allow the merging processes of each sample only wait for relative children processes from the previous step, instead of waiting for all the processes to be done
// change from .concat to .mix because .concat will wait for all the items proceeding from the first channel were emitted
(mergedChannelSomatic, mergedChannelGermline) = aMergedChannel.mix( iMergedChannel, wMergedChannel).map{
  item ->
    def key = item[2]+"__"+item[3]+"@"+item[0] // adding one unique key
    def target = item[0]
    def assay = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def tumorBam = item[4]
    def normalBam = item[5]
    def tumorBai = item[6]
    def normalBai = item[7]
    def intervalBed = item[8]

    return [ key, target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed ]
}.map{ 
    key, target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed -> 
    tuple ( 
         groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
         target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed
    )
}.transpose().into(2)


// if using strelka2, one should have manta for small InDels

if('strelka2' in tools) {
  tools.add('manta')
}

// --- Run Delly
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDelly, bamFiles) = bamFiles.into(2)

process SomaticDellyCall {
  tag {idTumor + "__" + idNormal + '@' + svType}

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}__${idNormal}_${svType}.filter.bcf") into dellyFilterOutput

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

// --- Run Mutect2
process RunMutect2 {
  tag {idTumor + "__" + idNormal + "@" + intervalBed.baseName}

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set id, target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelSomatic 
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
forMutect2Combine = forMutect2Combine.groupTuple()


// Combine Mutect2 VCFs, bcftools
process SomaticCombineMutect2Vcf {
  tag {idTumor + "__" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/mutations/mutect2", mode: params.publishDirMode }

  input:
    set id, idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcfOutput

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

(bamsForManta, bamsForStrelka, bamFiles) = bamFiles.into(3)

// --- Run Manta
process SomaticRunManta {
  tag {idTumor + "__" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/structural_variants/manta", mode: params.publishDirMode, pattern: "*.manta.vcf.{gz,gz.tbi}" }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.manta.vcf.gz") into mantaOutput
    set idTumor, idNormal, target, file("${outputPrefix}.manta.vcf.gz.tbi") into mantatbi
    set idTumor, idNormal, target, assay, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: "manta" in tools && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  options = ""
  if (assay == "wes") options = "--exome"
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

dellyFilterOutput = dellyFilterOutput.groupTuple(by: [0,1,2], size: 5)

dellyMantaCombineChannel = dellyFilterOutput.combine(mantaOutput, by: [0,1,2])

// --- Process Delly and Manta VCFs 
(sampleIdsForDellyMantaMerge, bamFiles) = bamFiles.into(2)

// Merge VCFs, Delly and Manta
process SomaticMergeDellyAndManta {
  tag {idTumor + "__" + idNormal}

  if (publishAll) {
    publishDir "${params.outDir}/somatic/structural_variants/delly", mode: params.publishDirMode, pattern: "*.delly.vcf.{gz,gz.tbi}"
  }

  input:
    set idTumor, idNormal, target, file(dellyBcfs), file(mantaFile) from dellyMantaCombineChannel

  output:
    file("${outputPrefix}.delly.manta.vcf.{gz,gz.tbi}") into vcfDellyMantaMergedOutput
    set file("${outputPrefix}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz"), file("${outputPrefix}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz.tbi") into somaticDellyVcfs

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

  if (publishAll) { publishDir "${params.outDir}/somatic/mutations/strelka2", mode: params.publishDirMode, pattern: "*.vcf.{gz,gz.tbi}"}

  input:
    set idTumor, idNormal, target, assay, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])

  output:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file('*strelka2.vcf.gz'), file('*strelka2.vcf.gz.tbi') into strelkaOutputMerged

  when: tools.containsAll(["manta", "strelka2"]) && runSomatic

  script:
  options = ""
  intervals = wgsTargets
  if (assay == "wes") {
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


mutectStrelkaChannel = mutect2CombinedVcfOutput.combine(strelkaOutputMerged, by: [0,1,2])

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
    set idTumor, idNormal, target, file("${outputPrefix}.pass.vcf") into vcfMergedOutput

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

  if (publishAll) {
    publishDir "${params.outDir}/somatic/mutations", mode: params.publishDirMode, pattern: "*.unfiltered.maf"
  }

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutput 
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFile
    file("${outputPrefix}.unfiltered.maf") into unfilteredMafFile

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

(bamsForMsiSensor, bamFiles) = bamFiles.into(2)

// --- Run MSIsensor
process RunMsiSensor {
  tag {idTumor + "__" + idNormal}

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.msisensor.tsv") into msiOutput

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

(msiOutputForMetaData, msiOutput) = msiOutput.into(2)

(bamFilesForSnpPileup, bamFiles) = bamFiles.into(2)

// --- Run FACETS 
process DoFacets {
  tag {idTumor + "__" + idNormal}

  // publishDir "${params.outDir}/somatic/facets", mode: params.publishDirMode, pattern: "*/*/*.Rdata"
  publishDir "${params.outDir}/somatic/facets/${tag}", mode: params.publishDirMode, pattern: "*.snp_pileup.dat.gz"
  publishDir "${params.outDir}/somatic/facets/${tag}", mode: params.publishDirMode, pattern: "${outputDir}/*.{Rdata,.png}"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamFilesForSnpPileup
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set assay, target, idTumor, idNormal, file("${outfile}") into SnpPileup
    set idTumor, idNormal, target, file("${outputDir}/*purity.out"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.seg"), file("${outputDir}/*hisens.out"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*purity.CNCF.png"), val("${outputFacetsSubdirectory}/${outputDir}") into FacetsOutput
    set file("${outputDir}/*purity.seg"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.CNCF.png"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.out") into FacetsPurity
    set file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.out") into FacetsHisens
    file("${tag}_OUT.txt") into FacetsPurityHisensOutput

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
    ${bamTumor} ${bamNormal}

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

(bamsForPolysolver, bamFiles) = bamFiles.into(2)

// Run Polysolver
process RunPolysolver {
  tag {idTumor + "__" + idNormal}
  
  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForPolysolver

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.hla.txt") into hlaOutput

  when: "polysolver" in tools && runSomatic
  
  script:
  outputPrefix = "${idTumor}__${idNormal}"
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


(bamsForLOHHLA, bamFiles) = bamFiles.into(2)

// Channel currently in order [ assay, target, tumorID, normalID, tumorBam, normalBam, tumorBai, normalBai ]

// Re-order bamsForLOHHLA into idTumor, idNormal, and target, i.e. 
// [ tumorID, normalID, target, tumorBam, normalBam, tumorBai, normalBai ]

bamsForLOHHLA = bamsForLOHHLA.map{ 
  item -> 
    def assay = item[0]
    def target = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def tumorBam = item[4]
    def normalBam = item[5]
    def tumorBai = item[6]
    def normalBai = item[7]

    return [ idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai ]
  }

// Polysolver channel currently in order []
// [ idTumor, idNormal, target, winners.hla.txt ]

// FACETS channel in order
// [ idTumor, idNormal, target, file("${outputDir}/*purity.Rdata"), file("${outputDir}/*.*") ]
// [idTumor, idNormal, target, *purity.out, *purity.cncf.txt, *purity.Rdata, purity.seg, hisens.out, hisens.cncf.txt, hisens.Rdata, hisens.seg into FacetsOutput


(facetsForLOHHLA, FacetsforMafAnno, FacetsOutput) = FacetsOutput.into(3)

facetsForLOHHLA = facetsForLOHHLA.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def purity_out = item[3]
    def purity_cncf = item[4]
    def purity_rdata = item[5]
    def purity_seg = item[6]
    def hisens_out = item[7]
    def hisens_cncf = item[8]
    def hisens_rdata = item[9]
    def hisens_seg = item[10]
    def purityCNCF_png = item[11]
    def hisensCNCF_png = item[12]

    return [ idTumor, idNormal, target, purity_out ]
  }

// create channel for MetaDataParser, using *purity.out from FACETS,

(facetsForLOHHLA, facetsForMetaDataParser) = facetsForLOHHLA.into(2)

(hlaOutputForLOHHLA, hlaOutput) = hlaOutput.into(2)

// *purity.out from FACETS, winners.hla.txt from POLYSOLVER, with the above

mergedChannelLOHHLA = bamsForLOHHLA.combine(hlaOutputForLOHHLA, by: [0,1,2]).combine(facetsForLOHHLA, by: [0,1,2])

// Run LOHHLA
process RunLOHHLA {
  tag {idTumor + "__" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/lohhla", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("winners.hla.txt"), file("*_purity.out") from mergedChannelLOHHLA
    set file(hlaFasta), file(hlaDat) from Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat])

  output:
    set file("*HLAlossPrediction_CI.txt"), file("*DNA.IntegerCPN_CI.txt"), file("*.pdf") optional true into lohhlaOutput

  when: tools.containsAll(["lohhla", "polysolver", "facets"]) && runSomatic

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  cat winners.hla.txt | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
  
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
    mv Figures/*.pdf .
  fi
  """
}

(mafFileForMafAnno, mafFileForMutSig, mafFile) = mafFile.into(3)

// --- Run Mutational Signatures, github.com/mskcc/mutation-signatures
process RunMutationSignatures {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.mutsig.txt") into mutSigOutput
    file("${outputPrefix}.mutsig.txt") into mutSigForAggregate

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

//Formatting the channel to be: idTumor, idNormal, target, purity_rdata

FacetsforMafAnno = FacetsforMafAnno.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def purity_out = item[3]
    def purity_cncf = item[4]
    def purity_rdata = item[5]
    def purity_seg = item[6]
    def hisens_out = item[7]
    def hisens_cncf = item[8]
    def hisens_rdata = item[9]
    def hisens_seg = item[10]
    def purityCNCF_png = item[11]
    def hisensCNCF_png = item[12]
    def facetsPath = item[13]
    
    return [idTumor, idNormal, target, purity_rdata, purity_cncf, hisens_cncf, facetsPath]
  }


(FacetsforMafAnno, FacetsforMafAnnoGermline) = FacetsforMafAnno.into(2)
facetsMafFileSomatic = FacetsforMafAnno.combine(mafFileForMafAnno, by: [0,1,2])


// --- Do FACETS MAF annotation and post processing
process SomaticFacetsAnnotation {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/facets/${facetsPath}", mode: params.publishDirMode, pattern: "*{armlevel,genelevel}.unfiltered.txt"
  
  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileSomatic

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.facets.maf"), file("${outputPrefix}.armlevel.unfiltered.txt") into FacetsAnnotationOutputs
    set file("${outputPrefix}.armlevel.unfiltered.txt"), file("${outputPrefix}.genelevel.unfiltered.txt") into FacetsArmGeneOutputs

  when: tools.containsAll(["facets", "mutect2", "manta", "strelka2"]) && runSomatic

  script:
  mapFile = "${idTumor}_${idNormal}.map"
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
  """
}

(mafFileForNeoantigen, FacetsAnnotationOutputs) = FacetsAnnotationOutputs.into(2)

//Formatting the channel to be: idTumor, idNormal, target, MAF

mafFileForNeoantigen = mafFileForNeoantigen.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def mafFile = item[3]
    def armLevel = item[4]
    return [idTumor, idNormal, target, mafFile]
  }

hlaOutput = hlaOutput.combine(mafFileForNeoantigen, by: [0,1,2])

(hlaOutputForMetaDataParser, hlaOutput) = hlaOutput.into(2)

// --- Run neoantigen prediction pipeline
process RunNeoantigen {
  tag {idTumor + "__" + idNormal}

  input:
    set idTumor, idNormal, target, file(polysolverFile), file(mafFile) from hlaOutput
    set file(neoantigenCDNA), file(neoantigenCDS) from Channel.value([
      referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS
    ])

  output:
    set idTumor, idNormal, target, file("${outputDir}/*") into neoantigenOut
    file("${idTumor}__${idNormal}.all_neoantigen_predictions.txt") into NetMhcStatsOutput
    file("${outputDir}/*.maf") into NeoantigenMafOutput

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
  outputDir = "neoantigen"
  tmpDir = "${outputDir}-tmp"
  tmpDirFullPath = "\$PWD/${tmpDir}/"  // must set full path to tmp directories for netMHC and netMHCpan to work; for some reason doesn't work with /scratch, so putting them in the process workspace
  """
  echo -e "${outputPrefix}\t`wc -l ${mafFile} | cut -d ' ' -f1`" > file-size.txt
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

// [idTumor, idNormal, target, armLevel]
FacetsAnnotationOutputs = FacetsAnnotationOutputs.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def mafFile = item[3]
    def armLevel = item[4]
    return [idTumor, idNormal, target, armLevel]
  }

mergedChannelMetaDataParser = facetsForMetaDataParser.combine(FacetsAnnotationOutputs, by: [0,1,2]).combine(msiOutputForMetaData, by: [0,1,2]).combine(hlaOutputForMetaDataParser, by: [0,1,2]).combine(mutSigOutput, by: [0,1,2]).unique()

// --- Generate sample-level metadata
process MetaDataParser {
  tag {idTumor + "__" + idNormal}
 
  input:
    set idTumor, idNormal, target, file(purityOut), file(armLevel), file(msifile), file(polysolverFile), file(mafFile), file(mutSigOutput) from mergedChannelMetaDataParser
    set file(idtCodingBed), file(agilentCodingBed), file(wgsCodingBed) from Channel.value([
      referenceMap.idtCodingBed, referenceMap.agilentCodingBed, referenceMap.wgsCodingBed
    ]) 

  output:
    file("*.sample_data.txt") into MetaDataOutputs

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
    --mutational_signatures_output ${mutSigOutput} \
    --polysolver_output ${polysolverFile} \
    --MAF_input ${mafFile} \
    --coding_baits_BED ${codingRegionsBed}
  
  mv ${idTumor}__${idNormal}_metadata.txt ${idTumor}__${idNormal}.sample_data.txt
  """
}

process SomaticAggregateMaf {
 
  publishDir "${params.outDir}/somatic", mode: params.publishDirMode

  input:
    file(mafFile) from NeoantigenMafOutput.collect()
    
  output:
    file("mut_somatic.maf") into MafFileOutput

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
 
  publishDir "${params.outDir}/somatic", mode: params.publishDirMode

  input:
    file(netmhcCombinedFile) from NetMhcStatsOutput.collect()

  output:
    file("mut_somatic_neoantigens.txt") into NetMhcChannel

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
 
  publishDir "${params.outDir}/somatic", mode: params.publishDirMode

  input:
    file(purityFiles) from FacetsPurity.collect()
    file(hisensFiles) from FacetsHisens.collect()
    file(purityHisensOutput) from FacetsPurityHisensOutput.collect()
    file(annotationFiles) from FacetsArmGeneOutputs.collect()

  output:
    set file("cna_hisens_run_segmentation.seg"), file("cna_purity_run_segmentation.seg") into FacetsMergedChannel
    set file("cna_armlevel.txt"), file("cna_genelevel.txt"), file("cna_facets_run_info.txt") into FacetsAnnotationMergedChannel
    
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
  cat facets_tmp/*armlevel.unfiltered.txt | grep -v "DIPLOID" | grep -v "Tumor_Sample_Barcode"  || [[ \$? == 1 ]] >> cna_armlevel.txt
  """
}

process SomaticAggregateSv {
 
  publishDir "${params.outDir}/somatic", mode: params.publishDirMode

  input:
    file(dellyMantaVcf) from vcfDellyMantaMergedOutput.collect()

  output:
    file("sv_somatic.vcf.{gz,gz.tbi}") into VcfBedPeChannel

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
 
  publishDir "${params.outDir}/somatic", mode: params.publishDirMode

  input:
    file(metaDataFile) from MetaDataOutputs.collect()

  output:
    file("sample_data.txt") into MetaDataOutputChannel

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

/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/

// GATK HaplotypeCaller
process GermlineRunHaplotypecaller {
  tag {idNormal + "@" + intervalBed.baseName}

  input:
    set id, target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set id, idTumor, idNormal, target, file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz.tbi"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz.tbi") into haplotypecallerOutput

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
haplotypecallerOutput = haplotypecallerOutput.groupTuple()

// merge VCFs, GATK HaplotypeCaller

process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/mutations/haplotypecaller", mode: params.publishDirMode }

  input:
    set id, idTumor, idNormal, target, file(haplotypecallerSnpVcf), file(haplotypecallerSnpVcfIndex), file(haplotypecallerIndelVcf), file(haplotypecallerIndelVcfIndex) from haplotypecallerOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools && runGermline 

  script: 
  idTumor = id.toString().split("__")[0]
  idNormal = id.toString().split("@")[0].split("__")[1]
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

(bamsForMantaGermline, bamsForStrelkaGermline, bamFiles) = bamFiles.into(3)

// --- Run Manta, germline
process GermlineRunManta {
  tag {idNormal}
  
  if (publishAll) { publishDir "${params.outDir}/germline/structural_variants/manta", mode: params.publishDirMode }
  
  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMantaGermline
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}.manta.vcf.gz"), file("${idNormal}.manta.vcf.gz.tbi") into mantaOutputGermline

  when: 'manta' in tools && runGermline

  // flag with --exome if exome
  script:
  options = ""
  if (assay == "wes") options = "--exome"
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

// --- Run Strelka2, germline
process GermlineRunStrelka2 {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/mutations/strelka2", mode: params.publishDirMode }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelkaGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsIntervals),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervalsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])
    
  output:
    set idTumor, idNormal, target, file("${idNormal}.strelka2.vcf.gz"), file("${idNormal}.strelka2.vcf.gz.tbi") into strelkaOutputGermline

  when: 'strelka2' in tools && runGermline
  
  script:
  options = ""
  intervals = wgsIntervals
  if (assay == "wes") {
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
haplotypecallerStrelkaChannel = haplotypecallerCombinedVcfOutput.combine(strelkaOutputGermline, by: [0,1,2])

(bamsForCombineChannel, bamFiles) = bamFiles.into(2)

bamsForCombineChannel = bamsForCombineChannel.map{
  item -> 
    def assay = item[0]
    def target = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def bamTumor = item[4]
    def bamNormal = item[5]
    def baiTumor = item[6]
    def baiNormal = item[7]
    
    return [idTumor, idNormal, target, assay, bamTumor, baiTumor]
  }

mergedChannelVcfCombine = bamsForCombineChannel.combine(haplotypecallerStrelkaChannel, by: [0,1,2])

// --- Combine VCFs with germline calls from Haplotypecaller and Strelka2
process GermlineCombineChannel {
  tag {idNormal}

  input:
    set idTumor, idNormal, target, assay, file(bamTumor), file(baiTumor), file(haplotypecallercombinedVcf), file(haplotypecallercombinedVcfIndex), file(strelkaVcf), file(strelkaVcfIndex) from mergedChannelVcfCombine
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
    set idTumor, idNormal, target, file("${idTumor}__${idNormal}.germline.vcf") into vcfMergedOutputGermline

  when: tools.containsAll(["strelka2", "haplotypecaller"]) && runGermline

  script:  
  isecDir = "${idNormal}.isec"
  gnomad = gnomadWgsVcf
  if (assay == 'wgs') {
    gnomad = gnomadWgsVcf
  }
  else if (assay == 'wes') {
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
  tag {idNormal}

  if (publishAll) {
    publishDir "${params.outDir}/germline/mutations", mode: params.publishDirMode, pattern: "*.unfiltered.maf"
  }

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutputGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFileGermline
    file("${outputPrefix}.unfiltered.maf") into unfilteredMafFileGermline

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

(mafFileGermline, mafFileGermlineFacets) = mafFileGermline.into(2)

facetsMafFileGermline = FacetsforMafAnnoGermline.combine(mafFileGermlineFacets, by: [0,1,2])

process GermlineFacetsAnnotation {
  tag {idNormal}

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), facetsPath, file(maf) from facetsMafFileGermline

  output:
    file("${outputPrefix}.facets.zygosity.maf") into mafFileAnnotatedGermline

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

  Rscript --no-init-file /usr/bin/annotate-with-zygosity-germline.R ${outputPrefix}.facets.maf ${outputPrefix}.facets.zygosity.maf
  """
}

// --- Call germline SVs with Delly
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDellyGermline, bamFiles) = bamFiles.into(2)

process GermlineDellyCall {
  tag {idNormal + '@' + svType}

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDellyGermline
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}_${svType}.filter.bcf") into dellyFilterOutputGermline

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

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)


dellyFilterOutputGermline = dellyFilterOutputGermline.groupTuple(by: [0,1,2], size: 5)

dellyMantaChannelGermline = dellyFilterOutputGermline.combine(mantaOutputGermline, by: [0,1,2])

// --- Merge Delly and Manta VCFs 
process GermlineMergeDellyAndManta {
  tag {idNormal}

  if (publishAll) {
    publishDir "${params.outDir}/germline/structural_variants/delly", mode: params.publishDirMode, pattern: "*delly.vcf.{gz,gz.tbi}"
  }

  input:
    set idTumor, idNormal, target, file(dellyBcf), file(mantaVcf), file(mantaVcfIndex) from dellyMantaChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}.delly.manta.vcf.gz"), file("${idNormal}.delly.manta.vcf.gz.tbi") into vcfFilterDellyMantaOutputGermline
    set file("${idNormal}.delly.manta.vcf.gz"), file("${idNormal}.delly.manta.vcf.gz.tbi") into germlineVcfBedPe
    set file("${idNormal}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz"), file("${idNormal}_{BND,DEL,DUP,INS,INV}.delly.vcf.gz.tbi") into germlineDellyVcfs

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


// --- Aggregate per-sample germline data, MAF
process GermlineAggregateMaf {

  publishDir "${params.outDir}/germline/", mode: params.publishDirMode

  input:
    file(mafFile) from mafFileAnnotatedGermline.collect()

  output:
    file("mut_germline.maf") into GermlineMafFileOutput
  
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

germlineVcfBedPe = germlineVcfBedPe.unique { new File(it.toString()).getName() }

// --- Aggregate per-sample germline data, SVs
process GermlineAggregateSv {
 
  publishDir "${params.outDir}/germline", mode: params.publishDirMode

  input:
    file(dellyMantaVcf) from germlineVcfBedPe.collect()

  output:
    file("sv_germline.vcf.{gz,gz.tbi}") into GermlineVcfBedPeChannel
  
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


/*
================================================================================
=                              Quality Control                                 =
================================================================================
*/


(bamsForPileup, bamFiles) = bamFiles.into(2)

allBamFiles = bamsForPileup.map{
  item ->
    def idTumor = item[2]
    def idNormal = item[3]
    def bamTumor = item[4]
    def bamNormal = item[5]
    def baiTumor = item[6]
    def baiNormal = item[7]

    return [[idTumor, idNormal], [bamTumor, bamNormal], [baiTumor, baiNormal]]
}.transpose().unique()

process QcPileup {
  tag {idSample}

  publishDir "${params.outDir}/qc/pileup/", mode: params.publishDirMode

  input:
    set idSample, file(bam), file(bai) from allBamFiles
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idSample, file("${idSample}.pileup") into (tumorPileup, normalPileup)

  when: !params.test && "pileup" in tools

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

(bamsForPileupTumor, bamsForPileupNormal, bamFiles) = bamFiles.into(3)

(tumorPileupConpairAll, tumorPileupConpair) = bamsForPileupTumor.combine(tumorPileup)
			 .filter{ item ->
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def bamNormal = item[5]
                            def baiTumor = item[6]
                            def baiNormal = item[7]
			    def idSample = item[8]
			    def samplePileup = item[9]
			    idSample == idTumor
                         }.map { item ->
                            def conpair = "conpair"
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def tumorPileup = item[9]

                            return [ conpair, idTumor, idNormal, tumorPileup ]
                         }.into(2)
(normalPileupConpairAll, normalPileupConpair) = bamsForPileupNormal.combine(normalPileup)
                         .filter{ item ->
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def bamNormal = item[5]
                            def baiTumor = item[6]
                            def baiNormal = item[7]
                            def idSample = item[8]
                            def samplePileup = item[9]
                            idSample == idNormal
                         }.map { item ->
                            def conpair = "conpair"
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def normalPileup = item[9]

                            return [ conpair, idTumor, idNormal, normalPileup ]
                         }.into(2)


pileupConpair = tumorPileupConpair.combine(normalPileupConpair, by: [0, 1, 2])

process QcConpair {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/qc/conpair/${idTumor}__${idNormal}", mode: params.publishDirMode

  input:
    set conpair, idTumor, idNormal, file(pileupTumor), file(pileupNormal) from pileupConpair
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    file("${outPrefix}.concordance.txt") into conpairConcordance
    file("${outPrefix}.contamination.txt") into conpairContamination

  when: !params.test && "conpair" in tools

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

pileupConpairAll = tumorPileupConpairAll.combine(normalPileupConpairAll, by: [0])

process QcConpairAll {
  tag {idTumor + "@" + idNormal}

  input:
    set conpair, idTumor, idNormal_noUse, file(pileupTumor), idTumor_noUse, idNormal, file(pileupNormal) from pileupConpairAll
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    file("${outPrefix}.concordance.txt") into conpairAllConcordance
    file("${outPrefix}.contamination.txt") into conpairAllContamination

  when: !params.test && "conpairall" in tools

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

process QcConpairAggregate {

  publishDir "${params.outDir}/qc/conpairAll/", mode: params.publishDirMode

  input:
    file(concordance) from conpairAllConcordance.collect()
    file(contamination) from conpairAllContamination.collect()

  output:
    set file('ConpairAll-concordance.txt'), file('ConpairAll-contamination.txt') into conpairAggregated

  when: !params.test

  script:
  """
  grep -v "concordance" *.concordance.txt | sed 's/.concordance.txt:/\t/' | cut -f1,3 | sort -k1,1 > ConpairAll-concordance.txt
  echo -e "Pairs\tSample_Type\tSample_ID\tContamination" > ConpairAll-contamination.txt
  grep -v "Contamination" *.contamination.txt | sed 's/.contamination.txt:/\t/' | sort -k1,1 >> ConpairAll-contamination.txt
  """
}



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



