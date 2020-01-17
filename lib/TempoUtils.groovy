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
      if(!checkDuplicatedRows(allRows, row, row, tsvFile)){System.exit(1)}

      [row.TUMOR_ID, row.NORMAL_ID]
    }
  }


  // No header need
  static def extractCohort(tsvFile) {
    def allRows = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t')
    .map { row ->
      if(!checkNumberOfItem(row, 1, tsvFile)){System.exit(1)}
      if(!checkDuplicatedRows(allRows, row, row, tsvFile)){System.exit(1)}
      def path = file(row[0], checkIfExists: true)
      def section = file(path.getParent()).getName().toString()

      [section, path]
    }
  }


  static def extractBAM(tsvFile, assayType) {
    def allRows = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkHeader([row.SAMPLE, row.TARGET, row.BAM, row.BAI], tsvFile)
      if(!checkNumberOfItem(row, 4, tsvFile)){System.exit(1)}
      if(!checkDuplicatedRows(allRows, row, row, tsvFile)){System.exit(1)}
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
    def allReadNames = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkHeader([row.SAMPLE, row.TARGET, row.FASTQ_PE1, row.FASTQ_PE2], tsvFile)
      if(!checkNumberOfItem(row, 4, tsvFile)){System.exit(1)}
      if(!checkDuplicatedRows(allRows, row, row, tsvFile)){System.exit(1)}
      if(!checkTarget(row.TARGET, assayType)){System.exit(1)}
      def idSample = row.SAMPLE
      def target = row.TARGET
      def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: true)
      def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: true)
      def fastqInfo = flowcellLaneFromFastq(fastqFile1)
      def fileID = fastqFile1.baseName.replaceAll("_+R1(?!.*R1)", "").replace(".fastq", "") + "@" + fastqInfo[0].replaceAll(":","@")
      
      if(!checkFileExtension(fastqFile1,".fastq.gz")){System.exit(1)}
      if(!checkFileExtension(fastqFile2,".fastq.gz")){System.exit(1)}

      def readName = fastqInfo[2]
      if(allReadNames.containsKey(readName)){
	println "ERROR: The following two files look like same file, because they contain the same read name ${readName} in the first line"
	println allReadNames.get(readName)
	println idSample + "\t" + row.FASTQ_PE1
	System.exit(1)
      }
      else{
        allReadNames[readName] = idSample + "\t" + row.FASTQ_PE1
      }

      [idSample, target, fileID, fastqFile1, fastqFile2]
    }
  }


 static def flowcellLaneFromFastq(path) {
    // https://github.com/SciLifeLab/Sarek/blob/917a4d7f4dceb5a524eb7bd1c287cd197febe9c0/main.nf#L639-L666
    // parse first line of a FASTQ file (optionally gzip-compressed)
    // and return the flowcell id and rgID number.
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    String fullName
    if (fields.size() == 7) {
      // CASAVA 1.8+ format
      // we include instrument name and run id in fcid to ensure the uniqueness
      fcid = fields[0] + ":" + fields[1] + ":" + fields[2]
      lane = fields[3].toInteger()
      fullName = "@" + line
    }
    else if (fields.size() == 5) {
      fcid = fields[0]
      lane = fields[1].toInteger()
      fullName = "@" + line
    }
    [fcid, lane, fullName]
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
    else if(assayType == "exome"){ supportedTargets = ["agilent", "idt"]}
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
  static def checkDuplicatedRows(hash, key, value, tsv) {
    if(hash.containsKey(key)){
        println "ERROR: Duplicatd inputs found in ${tsv}:"
        println hash.get(key)
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
