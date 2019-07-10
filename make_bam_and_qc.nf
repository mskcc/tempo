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

if (params.mapping) mappingPath = params.mapping
if (params.pairing) pairingPath = params.pairing

if (!check_for_duplicated_rows(pairingPath)) {
  println "ERROR: Duplicated row found in pairing file. Please fix the error and rerun the pipeline"
  exit 1
}

if (!check_for_mixed_assay(mappingPath)) {
  println "ERROR: You can only use either assays 'exome' or 'genome', not both WES and WGS together"
  exit 1
}

outname = params.outname

referenceMap = defineReferenceMap()

fastqFiles = Channel.empty()

mappingFile = file(mappingPath)
pairingfile = file(pairingPath)

pairingT = extractPairing(pairingfile)

fastqFiles = extractFastq(mappingFile)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


fastqFiles.groupTuple(by:[0]).map { key, lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets -> tuple( groupKey(key, lanes.size()), lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets) }.set { groupedFastqs }

groupedFastqs.into { groupedFastqsDebug; fastPFiles; fastqFiles }
fastPFiles = fastPFiles.transpose()
fastqFiles = fastqFiles.transpose()

// AlignReads - Map reads with BWA mem output SAM

process AlignReads {
  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  input:
    set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assay, targetFile from fastqFiles
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idSample, lane, file("${lane}.bam"), assay, targetFile into (unsortedBam)

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
    set idSample, lane, file("${lane}.bam"), assay, targetFile from unsortedBam

  output:
    set idSample, lane, file("${lane}.sorted.bam"), assay, targetFile into (sortedBam, sortedBamDebug)

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

sortedBam.groupTuple().set { groupedBam }
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
groupedBam.into { groupedBamDebug; groupedBam }

process MergeBams {
  tag {idSample}

  input:
    set idSample, lane, file(bam), assay, targetFile from groupedBam

  output:
    set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile into (mergedBam, mergedBamDebug)

  script:
  """
  samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
  """
}

// GATK MarkDuplicates

process MarkDuplicates {
  tag {idSample}

   publishDir "${params.outDir}/MarkDup/${idSample}", mode: params.publishDirMode

  input:
    set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile from mergedBam

  output:
    set file("${idSample}.md.bam"), file("${idSample}.md.bai"), idSample, lane, assay, targetFile into duplicateMarkedBams
    set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), assay, targetFile into markDuplicatesTSV
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
    bam, bai, idSample, lane, assay, targetFile ->
    tag = bam.baseName.tokenize('.')[0]
    [idSample, bam, bai, assay, targetFile]
}

(mdBam, mdBamToJoin, mdDebug) = duplicateMarkedBams.into(3)

process CreateRecalibrationTable {
  tag {idSample}

  input:
    set idSample, file(bam), file(bai), assay, targetFile from mdBam 

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
    set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), val("${idSample}.recal.table"), assay, targetFile into recalibrationTableTSV

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

recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0])

// GATK ApplyBQSR, RecalibrateBAM

process RecalibrateBam {
  tag {idSample}

  publishDir "${params.outDir}/BQSR/${idSample}", mode: params.publishDirMode

  input:
    set idSample, file(bam), file(bai), assay, targetFile, file(recalibrationReport) from recalibrationTable
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict 
    ])

  output:
    set idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai"), assay, targetFile into recalibratedBam, recalibratedBamForStats, recalibratedBamForOutput, recalibratedBamForOutput2
    set idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bam.bai"), assay, targetFile into recalibratedBamTSV
    val(idSample) into currentSample
    file("${idSample}.recal.bam") into currentBam
    file("${idSample}.recal.bai") into currentBai
    val(assay) into assays
    val(targetFile) into targets

  script:
  """
  gatk ApplyBQSR \
    --reference ${genomeFile} \
    --create-output-bam-index true \
    --bqsr-recal-file ${recalibrationReport} \
    --input ${bam} \
    --output ${idSample}.recal.bam

  cp -p ${idSample}.recal.bai ${idSample}.recal.bam.bai
  """
}


recalibratedBamForOutput.combine(pairingT)
                        .filter { item -> // only keep combinations where sample is same as tumor pair sample
                          def sampleID = item[0]
                          def sampleBam = item[1]
                          def sampleBai = item[2]
                          def assay = item[3]
                          def target = item[4]
                          def tumorID = item[5]
                          def normalID = item[6]
                          sampleID == tumorID
                        }.map { item -> // re-order the elements
                          def sampleID = item[0]
                          def sampleBam = item[1]
                          def sampleBai = item[2]
                          def assay = item[3]
                          def target = item[4]
                          def tumorID = item[5]
                          def normalID = item[6]
                          def tumorBam = sampleBam
                          def tumorBai = sampleBai

                          return [ assay, target, tumorID, normalID, tumorBam, tumorBai ]
                        }.combine(recalibratedBamForOutput2)
                        .filter { item ->
                          def assay = item[0]
                          def target = item[1]
                          def tumorID = item[2]
                          def normalID = item[3]
                          def tumorBam = item[4]
                          def tumorBai = item[5]
                          def sampleID = item[6]
                          def normalBam = item[7]
                          def normalBai = item[8]
                          sampleID == normalID
                        }.map { item -> // re-order the elements
                          def assay = item[0]
                          def target = item[1]
                          def tumorID = item[2]
                          def normalID = item[3]
                          def tumorBam = item[4]
                          def tumorBai = item[5]
                          def sampleID = item[6]
                          def normalBam = item[7]
                          def normalBai = item[8]

                          return [ assay, target, tumorID, normalID, tumorBam, normalBam, tumorBai, normalBai ]
                        }
                        .set { result }

File file = new File(outname)
file.newWriter().withWriter { w ->
    w << "ASSAY\tTARGET\tTUMOR_ID\tNORMAL_ID\tTUMOR_BAM\tNORMAL_BAM\tTUMOR_BAI\tNORMAL_BAI\n"
}

if (workflow.profile == 'awsbatch') {
    result.subscribe { Object obj ->
      file.withWriterAppend { out ->
        out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\ts3:/${obj[4]}\ts3:/${obj[5]}\ts3:/${obj[6]}\ts3:/${obj[7]}"
    }
  }
}
  else {
    result.subscribe { Object obj ->
      file.withWriterAppend { out ->
        out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\t${obj[4]}\t${obj[5]}\t${obj[6]}\t${obj[7]}"
    }
  }
}

// FastP - FastP on lane pairs, R1/R2

process FastP {
  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  publishDir "${params.outDir}/FastP/${idSample}", mode: params.publishDirMode

  input:
    set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assays, targetFiles from fastPFiles
    
  output:
    file("*.html") into fastPResults 

  script:
  """
  fastp -h ${lane}.html -i ${fastqFile1} -I ${fastqFile2}
  """
}


ignore_read_groups = Channel.from( true , false )

process Alfred {
  tag {idSample}

  publishDir "${params.outDir}/Alfred/${idSample}", mode: params.publishDirMode
  
  input:
    each ignore_rg from ignore_read_groups
    set idSample, file(bam), file(bai), assay, targetFile from recalibratedBam

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
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex")
  ]
}

def debug(channel) {
  channel.subscribe { Object obj ->
    println "DEBUG: ${obj.toString()};"
  }
}

def extractPairing(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    [row.NORMAL_ID, row.TUMOR_ID]
  }
}

def extractFastq(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    checkNumberOfItem(row, 6)
    def idSample = row.SAMPLE
    def lane = row.LANE
    def assayValue = row.ASSAY
    def targetFile = row.TARGET
    def fastqFile1 = returnFile(row.FASTQ_PE1)
    def sizeFastqFile1 = fastqFile1.size()
    def fastqFile2 = returnFile(row.FASTQ_PE2)
    def sizeFastqFile2 = fastqFile2.size()

    def assay = assayValue.toLowerCase() //standardize genome/wgs/WGS to wgs, exome/wes/WES to wes

    if ((assay == "genome") || (assay == "wgs")) {
      assay = "wgs"
    }
    if ((assay == "exome") || (assay == "wes")) {
      assay = "wes"
    }

    checkFileExtension(fastqFile1,".fastq.gz")
    checkFileExtension(fastqFile2,".fastq.gz")

    [idSample, lane, fastqFile1, sizeFastqFile1, fastqFile2, sizeFastqFile2, assay, targetFile]
  }
}

// Check file extension
def checkFileExtension(it, extension) {
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
  if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
  return true
}

// Return file if it exists
def returnFile(it) {
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  return file(it)
}

def check_for_duplicated_rows(pairingFilePath) {
  def entries = []
  file( pairingFilePath ).eachLine { line ->
    entries << line
  }
  return entries.toSet().size() == entries.size()
}

def check_for_mixed_assay(mappingFilePath) {
  def wgs = false
  def wes = false
  file( mappingFilePath ).eachLine { line ->
    currentLine = line.toLowerCase()
    if (currentLine.contains('\tgenome\t') || currentLine.contains('\twgs\t')) {
      wgs = true
    }
    if (currentLine.contains('\texome\t') || currentLine.contains('\twes\t')) {
      wes = true
    }
  return !(wgs && wes)
  }
}
