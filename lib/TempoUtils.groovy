import static nextflow.Nextflow.file
import nextflow.Channel

// NOTE: 
// static methods in Groovy are meant to be accessed directly from the class
// e.g. TempoUtils.debug()
//

class TempoUtils {

  static def extractPairing(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      [row.TUMOR_ID, row.NORMAL_ID]
    }
  }

  static def extractCohort(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      def path = file(row.PATH, checkIfExists: true)
      def section = file(path.getParent()).getName().toString()
      [section, path]
    }
  }


  static def extractBAM(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
//      checkNumberOfItem(row, 5)	// Disable check columns for now to support older version of input files, especially for Travis-CI
      def idSample = row.SAMPLE
      def target = row.TARGET
      def bam = file(row.BAM, checkIfExists: true)
      // check if using bam.bai or bam.bam.bai
      def bai = file(validateBamIndexFormat(row.BAM), checkIfExists: true)

      [idSample, target, file(bam), file(bai)]
    }
  }


  static def extractFastq(tsvFile) {
    def allReadNames = [:]
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
//      checkNumberOfItem(row, 4)	// Disable check columns for now to support older version of input files, especially for Travis-CI
      def idSample = row.SAMPLE
      def target = row.TARGET
      def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: true)
      def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: true)
      def fastqInfo = flowcellLaneFromFastq(fastqFile1)
      def fileID = fastqFile1.baseName.replaceAll("_+R1(?!.*R1)", "").replace(".fastq", "") + "@" + fastqInfo[0].replaceAll(":","@")
      
      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

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

      [idSample, fileID, fastqFile1, fastqFile2, target]
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

  // Check which format of BAM index used, input 'it' as BAM file 'bamTumor.bam'
  // not a static method, as currently written
  static def validateBamIndexFormat(it) {
    def bamFilename = it.take(it.lastIndexOf('.'))
    // Check BAM index extension
    if (file(bamFilename + ".bai").exists()){
      return("${bamFilename}.bai")
    } 
    else if (file(bamFilename + ".bam.bai").exists()){
      return("${bamFilename}.bam.bai")
    } 
    else {
      println "ERROR: Cannot find BAM indices for ${it}. Please index BAMs in the same directory with 'samtools index' and re-run the pipeline."
      System.exit(1)
    }
  }

  
  static def checkTargetAndAssayType(filePath, assayType) {
    def supportedTargets = []
    if(assayType == "genome"){ supportedTargets = ["wgs"] }
    else if(assayType == "exome"){ supportedTargets = ["agilent", "idt"]}
    else {println "ERROR: Unsupported --assayType ${assayType}. Supported values are \"exome\" and \"genome\"";System.exit(1)}

    def targetList = Channel.from(filePath).splitCsv(sep: '\t', header: true).map{row -> [row.TARGET]}.flatten().unique().toSortedList().get()
    def unsupported = targetList - supportedTargets

    if(unsupported != []){
      println "ERROR: Following targets are not supported by --assayType ${assayType}, only ${supportedTargets} supported. Please check your mapping file."
      println "${unsupported}"
      System.exit(1)
    }
  }
  
  
  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) {
	println "File: ${it} has the wrong extension: ${extension} see --help for more information"
	System.exit(1)
    }
  }

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number){
	println "Malformed row in TSV file: ${row}, see --help for more information"
	System.exit(1)
    }
    return true
  }

  static def check_for_duplicated_rows(FilePath) {
    def entries = []
    file( FilePath, checkIfExists: true ).eachLine { line ->
      if (!line.isEmpty()){
        entries << line
      }
    }
    if(entries.toSet().size() != entries.size()){
      println "ERROR: Duplicated row found in ${FilePath}. Please fix the error and re-run the pipeline."
      System.exit(1)
    }
  }

  static def crossValidateSamples(mapping, pairing) {
    def samplesInMapping = Channel.from(mapping).splitCsv(sep: '\t', header: true).map{[it.SAMPLE]}.flatten().unique().toSortedList().get()
    def samplesInPairing = extractPairing(pairing).flatten().unique().toSortedList().get()
    def mappingOnly = samplesInMapping - samplesInPairing
    def pairingOnly = samplesInPairing - samplesInMapping
    def extraSamples = mappingOnly + pairingOnly

    if(extraSamples != []){
	println "ERROR: The following samples are present in either mapping file or pairing file only. Please ensure all samples are present in both files."
        println "Mapping: ${mappingOnly}"
        println "Pairing: ${pairingOnly}"
	System.exit(1)
    }
  }

  static def crossValidateTargets(mappingFile, pairingFile) {
    Channel.from(mappingFile)
      .splitCsv(sep: '\t', header: true)
      .map { row ->
        def idSample = row.SAMPLE
	def target = row.TARGET
	[idSample, target]
      }
      .unique()
      .combine(extractPairing(pairingFile).unique())
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
