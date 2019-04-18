#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - AlignReads - Map reads with BWA mem output SAM
 - SortBAM - Sort BAM with samtools
 - MergeBam - Merge BAM for the same samples from different lanes
 - MarkDuplicates - Mark Duplicates with GATK4
 - CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
 - RecalibrateBam - Recalibrate Bam with PrintReads
*/

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

referenceMap = defineReferenceMap()

fastqFiles = Channel.empty()

tsvFile = file(tsvPath)

fastqFiles = (tsvFile) 

// Duplicate channel
fastqFiles.into { fastqFiles; fastQCFiles; fastPFiles }

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// FastP - FastP on lane pairs, R1/R2

process FastP {
  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  publishDir params.outDir, mode: params.publishDirMode

  input:
    set idSample, lane, file(fastqFile1), file(fastqFile2) from fastPFiles

  output:
    file("*.html") into fastPResults 

  """
  fastp -h ${lane}.html -i ${fastqFile1} -I ${fastqFile2}
  """
}

// AlignReads - Map reads with BWA mem output SAM

process AlignReads {
  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  input:
    set idSample, lane, file(fastqFile1), file(fastqFile2) from fastqFiles
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idSample, lane, file("${lane}.bam") into (unsortedBam)

  script:
    readGroup = "@RG\\tID:${lane}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
    
  """
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} | samtools view -Sb - > ${lane}.bam
  """
}

// SortBAM - Sort unsorted BAM with samtools, 'samtools sort'

process SortBAM {
  tag {lane}

  input:
    set idSample, lane, file("${lane}.bam") from unsortedBam

  output:
    set idSample, lane, file("${lane}.sorted.bam") into (sortedBam, sortedBamDebug)

  script:
  // Refactor when https://github.com/nextflow-io/nextflow/pull/1035 is merged
  if(params.mem_per_core) { 
    mem = task.memory.toString().split(" ")[0].toInteger() - 1 
  }
  else {
    mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
  }
  """
  samtools sort -m ${mem}G -@ ${task.cpus} -o ${lane}.sorted.bam ${lane}.bam
  """
}

singleBam = Channel.create()
singleBamDebug = Channel.create()
groupedBam = Channel.create()
groupedBamDebug = Channel.create()
sortedBam.groupTuple(by:[0,1])
  .choice(singleBam, groupedBam) {it[1].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
  idSample, lane, bam ->
  [idSample, bam]
}
sortedBamDebug.groupTuple(by:[0,1])
  .choice(singleBamDebug, groupedBamDebug) {it[1].size() > 1 ? 1 : 0}
singleBamDebug = singleBamDebug.map {
  idSample, lane, bam ->
  [idSample, bam]
}

if (params.debug) {
  debug(groupedBamDebug);
  debug(singleBamDebug);
}   

process MergeBams {
  tag {idSample}

  input:
    set idSample, lane, file(bam) from groupedBam

  output:
    set idSample, lane, file("${idSample}.merged.bam") into (mergedBam, mergedBamDebug)

  script:
  """
  samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
  """
}

if (params.debug) {
  debug(mergedBamDebug);
}

if (params.verbose) singleBam = singleBam.view {
  "Single BAM:\n\
  ID    : Sample: ${it[0]}\tLane: ${it[1]}\t\n\
  File  : [${it[2].fileName}]"
}

if (params.verbose) mergedBam = mergedBam.view {
  "Merged BAM:\n\
  ID    : Sample: ${it[0]}\tLane: ${it[1]}\t\n\
  File  : [${it[2]}]"
}

mergedBam = mergedBam.mix(singleBam)

if (params.verbose) mergedBam = mergedBam.view {
  "BAM for MarkDuplicates:\n\
  ID    : Sample: ${it[0]}\tLane: ${it[1]}\t\n\
  File  : [${it[2]}]"
}

// GATK MarkDuplicates

process MarkDuplicates {
  tag {idSample}

   publishDir params.outDir, mode: params.publishDirMode

  input:
    set idSample, lane, file("${idSample}.merged.bam") from mergedBam

  output:
    set file("${idSample}.md.bam"), file("${idSample}.md.bai"), idSample, lane into duplicateMarkedBams
    set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai") into markDuplicatesTSV
    file ("${idSample}.bam.metrics") into markDuplicatesReport

  script:
  """
  gatk MarkDuplicates --java-options ${params.markdup_java_options}  \
    --MAX_RECORDS_IN_RAM 50000 \
    --INPUT ${idSample}.merged.bam \
    --METRICS_FILE ${idSample}.bam.metrics \
    --TMP_DIR . \
    --ASSUME_SORT_ORDER coordinate \
    --CREATE_INDEX true \
    --OUTPUT ${idSample}.md.bam
  """
}

duplicateMarkedBams = duplicateMarkedBams.map {
    bam, bai, idSample, lane ->
    tag = bam.baseName.tokenize('.')[0]
    [idSample, bam, bai]
}

(mdBam, mdBamToJoin) = duplicateMarkedBams.into(2)
/*
if (params.verbose) mdBamToJoin = mdBamToJoin.view {
  "MD Bam to Join BAM:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}

if (params.verbose) mdBam = mdBam.view {
  "BAM for MarkDuplicates:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[3].fileName}]"
}
*/
process CreateRecalibrationTable {
  tag {idSample}

  input:
    set idSample, file(bam), file(bai) from mdBam 

    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex)  from Channel.value([
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
    set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), val("${idSample}.recal.table") into recalibrationTableTSV

  script:
  known = knownIndels.collect{ "--known-sites ${it}" }.join(' ')

  """
  gatk BaseRecalibrator \
    --tmp-dir /tmp \
    --reference ${genomeFile} \
    --known-sites ${dbsnp} \
    ${known} \
    --verbosity INFO \
    --input ${bam} \
    --output ${idSample}.recal.table
  """
}

recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0, 1])

process RecalibrateBam {
  tag {idSample}

  publishDir params.outDir, mode: params.publishDirMode

  input:
    set idSample, file(bam), file(bai), file(recalibrationReport) from recalibrationTable

    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict 
    ])

  output:
    set idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai") into recalibratedBam, recalibratedBamForStats
    set idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai") into recalibratedBamTSV

  script:
  """
  gatk ApplyBQSR \
    --reference ${genomeFile} \
    --create-output-bam-index true \
    --bqsr-recal-file ${recalibrationReport} \
    --input ${bam} \
    --output ${idSample}.recal.bam
  """
}

ignore_read_groups = Channel.from( true , false )

process Alfred {
  tag {idSample}

  publishDir params.outDir, mode: params.publishDirMode
  
  input:
    each ignore_rg from ignore_read_groups
    set idSample, file(bam), file(bai) from recalibratedBam

    file(genomeFile) from Channel.value([
      referenceMap.genomeFile
    ])

  output:
    set ignore_rg, idSample, file("*.tsv.gz"), file("*.tsv.gz.pdf") into bamsQCStats

  script:
  def ignore = ignore_rg ? "--ignore" : ''
  def outfile = ignore_rg ? "${idSample}.alfred.tsv.gz" : "${idSample}.alfred.RG.tsv.gz"
  """
  alfred qc --reference ${genomeFile} ${ignore} --outfile ${outfile} ${bam} && Rscript /opt/alfred/scripts/stats.R ${outfile}
  """

}

/*
================================================================================
=                               AWESOME FUNCTIONS                             =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  return [
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex'         : checkParamReturnFile("bwaIndex"), 
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
  ]
}

def debug(channel) {
  channel.subscribe { Object obj ->
    println "DEBUG: ${obj.toString()};"
  }
}

def extractFastq(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 4)
    def idSample =    row[0]
    def lane =   row[1]
    def fastqFile1 = SarekUtils.returnFile(row[2])
    def sizeFastqFile1 = fastqFile1.size()
    def fastqFile2 = SarekUtils.returnFile(row[3])
    def sizeFastqFile2 = fastqFile2.size()

    SarekUtils.checkFileExtension(fastqFile1,".fastq.gz")
    SarekUtils.checkFileExtension(fastqFile2,".fastq.gz")

    [idSample, lane, fastqFile1, sizeFastqFile1, fastqFile2, sizeFastqFile2]
    // [idSample, lane, fastqFile1, fastqFile2]
  }
}
