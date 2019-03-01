#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - AlignReads - Map reads with BWA mem output SAM
 - ConvertSAMtoBAM - Convert SAM to BAM with samtools
 - SortBAM - Sort BAM with samtools
 - MarkDuplicates - Mark Duplicates with GATK4
 - CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
 - RecalibrateBam - Recalibrate Bam with PrintReads
*/


/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/




// tsvPath 


// Evan comment: this is a silly thing to keep, 
// but it reminds me that this could be so much more flexible with complex conditionals

tsvPath = ''
if (params.sample) tsvPath = params.sample

referenceMap = defineReferenceMap()

fastqFiles = Channel.empty()

tsvFile = file(tsvPath)

fastqFiles = extractFastq(tsvFile)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


// tag
// https://www.nextflow.io/docs/latest/process.html#tag
// The tag directive allows you to associate each process executions with a custom label, 
// so that it will be easier to identify them in the log file or in the trace execution report.



// AlignReads - Map reads with BWA mem output SAM

process AlignReads {
  tag {idPatient + "-" + idRun}   // The tag directive allows you to associate each process executions with a custom label

  input:
    set idPatient, status, idSample, idRun, file(fastqFile1), file(fastqFile2) from fastqFiles
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idPatient, status, idSample, idRun, file("${idRun}.sam") into (alignedSam)

  script:
  readGroup = "@RG\\tID:Seq01p\\tSM:Seq01\\tPL:ILLUMINA\\tPI:330"
  """
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} > ${idRun}.sam
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
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
  ]
}


def extractFastq(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "subject gender status sample lane fastq1 fastq2"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 7)
    def idPatient  = row[0]
    def gender     = row[1]
    def status     = SarekUtils.returnStatus(row[2].toInteger())
    def idSample   = row[3]
    def idRun      = row[4]
    def fastqFile1 = SarekUtils.returnFile(row[5])
    def fastqFile2 = SarekUtils.returnFile(row[6])

    SarekUtils.checkFileExtension(fastqFile1,".fastq.gz")
    SarekUtils.checkFileExtension(fastqFile2,".fastq.gz")

    [idPatient, gender, status, idSample, idRun, fastqFile1, fastqFile2]
  }
}


