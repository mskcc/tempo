import static nextflow.Nextflow.file
import nextflow.Channel

// NOTE: 
// static methods in Groovy are meant to be accessed directly from the class
//

class TempoUtils {

  static def extractPairing(tsvFile) {
    def allRows = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkHeader([row.TUMOR_ID, row.NORMAL_ID], tsvFile)
      if(!checkNumberOfItem(row, 2, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allRows, row, row, tsvFile)){System.exit(1)}

      [row.TUMOR_ID, row.NORMAL_ID]
    }
  }


  // No header need
  static def extractCohort(tsvFile) {
    def allRows = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      if(!checkNumberOfItem(row, 3, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allRows, row, row, tsvFile)){System.exit(1)}
      checkHeader([row.TUMOR_ID, row.NORMAL_ID, row.COHORT], tsvFile)

      [row.COHORT, row.TUMOR_ID, row.NORMAL_ID, row.PATH]
    }
  }


  static def extractBAM(tsvFile, assayType) {
    def allRows = [:]
    def allFiles = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkHeader([row.SAMPLE, row.TARGET, row.BAM, row.BAI], tsvFile)
      if(!checkNumberOfItem(row, 4, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allRows, row, row, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allFiles, row.BAM, row.SAMPLE + "\t" + row.BAM, tsvFile)){System.exit(1)}
      if(!checkTarget(row.TARGET, assayType)){System.exit(1)}
      def idSample = row.SAMPLE
      def target = row.TARGET
      def bam = file(row.BAM, checkIfExists: true)
      if(!checkFileExtension(bam,".bam")){System.exit(1)}
      def bai = file(row.BAI, checkIfExists: true)
      if(!checkFileExtension(bai,".bai")){System.exit(1)}

      [idSample, target, file(bam), file(bai)]
    }
  }


  static def extractFastq(tsvFile, assayType) {
    def allRows = [:]
    def allFiles = [:]

    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkHeader([row.SAMPLE, row.TARGET, row.FASTQ_PE1, row.FASTQ_PE2], tsvFile)
      if(!checkNumberOfItem(row, 4, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allRows, row, row, tsvFile)){System.exit(1)}
      if(!checkDuplicates(allFiles, row.FASTQ_PE1, row.SAMPLE + "\t" + row.FASTQ_PE1, tsvFile)){System.exit(1)}
      if(!checkTarget(row.TARGET, assayType)){System.exit(1)}
      def idSample = row.SAMPLE
      def target = row.TARGET
      def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: true)
      def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: true)
      
      if(!checkFileExtension(fastqFile1,".fastq.gz")){System.exit(1)}
      if(!checkFileExtension(fastqFile2,".fastq.gz")){System.exit(1)}

      [idSample, target, fastqFile1, fastqFile2]
    }
  }


  // Check header
  static def checkHeader(row, tsv) {
    if(row.contains(null)){
      println "ERROR: Wrong header in ${tsv}. See manual for more infomation."
      System.exit(1)
    }
  }


  // Check supported assayType ("wgs" or "wes" are not supported, use "genome" or "exome" instead)
  static def checkAssayType(assayType) {
    if(assayType != "genome" && assayType != "exome"){
      println "ERROR: Unsupported \"--assayType ${assayType}\". Supported values are \"exome\" and \"genome\""
      System.exit(1)
    }
  }

  // Check TARGET field in mapping file is in the supported bait set list associated with --assayType what was provided as parameter
  static def checkTarget(it, assayType) {
    def supportedTargets = []
    if(assayType == "genome"){ supportedTargets = ["wgs"] }
    else if(assayType == "exome"){ supportedTargets = ["agilent", "idt","idt_v2"]}
    else {} // this is covered by checkAssayType(){} above

    if(!supportedTargets.contains(it)){
      println "ERROR: \"${it}\" is not a supported target (only ${supportedTargets}) for \"--assayType ${assayType}\". Please check your mapping file."
      return false
    }
    else{ return true }
  }
  
  
  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) {
	println "File: ${it} has the wrong extension: ${extension}. See manual for more information"
        return false
    }
    else{ return true }
  }

  // Check if a row has the minimal number of item
  static def checkNumberOfItem(row, number, tsv) {
    if (row.size() < number){
	println "Missing field (null) in the following row from ${tsv}: ${row}"
	return false
    }
    else{ return true }
  }

  // Check duplicate rows
  static def checkDuplicates(hash, key, value, tsv) {
    if(hash.containsKey(key)){
        println "ERROR: Duplicatd inputs found in ${tsv}"
	println ""
        println hash.get(key)
	println "${value}"
	return false
    }
    else{
        hash[key] = value
	return true
    }
  }


  // Check samples are present both mapping and pairing files
  static def crossValidateSamples(mapping, pairing) {
    def samplesInMapping = mapping.map{[it[0]]}.flatten().unique().toSortedList().get()
    def samplesInPairing = pairing.flatten().unique().toSortedList().get()
    def mappingOnly = samplesInMapping - samplesInPairing
    def pairingOnly = samplesInPairing - samplesInMapping
    def extraSamples = mappingOnly + pairingOnly

    if(extraSamples != []){
	println "ERROR: The following samples are present in either mapping file or pairing file only. Please ensure all samples are present in both files."
        println "Mapping: ${mappingOnly}"
        println "Pairing: ${pairingOnly}"
        return false
    }
    else{ return true }
  }


  // Check Tumor and Normal samples paired together have the same bait set
  static def crossValidateTargets(mapping, pairing) {
    mapping.map{[it[0], it[1]]}
      .unique()
      .combine(pairing)
      .filter{ item ->
	def mappingID = item[0]
	def target = item[1]
	def tumorID = item[2]
	def normalID = item[3]
	(mappingID == tumorID) || (mappingID == normalID)
      }
      .groupTuple(by: [2,3])
      .map{ item ->
        if(item[1].unique().size() > 1){
	  println "ERROR: Tumor and Normal pair ${item[0]} used differnt bait sets: ${item[1]}"
	  System.exit(1)
	}
      }
  }

}
