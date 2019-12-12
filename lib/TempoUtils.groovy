import static nextflow.Nextflow.file
import nextflow.Channel

// NOTE: 
// static methods in Groovy are meant to be accessed directly from the class
// e.g. TempoUtils.debug()
//

class TempoUtils {

  static def debug(channel) {
    channel.subscribe { Object obj ->
      println "DEBUG: ${obj.toString()};"
    }
  }

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
      [row.ID, row.PATH]
    }
  }

  static def extractFastq(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
//      checkNumberOfItem(row, 4)	// Disable check columns for now to support older version of input files, especially for Travis-CI
      def idSample = row.SAMPLE
      def targetFile = row.TARGET
      def fastqFile1 = returnFile(row.FASTQ_PE1)
      def fastqFile2 = returnFile(row.FASTQ_PE2)
      def fileID = fastqFile1.baseName.replaceAll("_+R1(?!.*R1)", "").replace(".fastq", "") + "@" + flowcellLaneFromFastq(fastqFile1)[0].replaceAll(":","@")
      
      checkFileExtension(fastqFile1,".fastq.gz")
      checkFileExtension(fastqFile2,".fastq.gz")

      [idSample, fileID, fastqFile1, fastqFile2, targetFile]
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
    if (fields.size() == 7) {
      // CASAVA 1.8+ format
      // we include instrument name and run id in fcid to ensure the uniqueness
      fcid = fields[0] + ":" + fields[1] + ":" + fields[2]
      lane = fields[3].toInteger()
    }
    else if (fields.size() == 5) {
      fcid = fields[0]
      lane = fields[1].toInteger()
    }
    [fcid, lane]
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
      exit 1
    }
  }

  
  static def extractBAM(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
//      checkNumberOfItem(row, 5)	// Disable check columns for now to support older version of input files, especially for Travis-CI
      def idTumor = row.TUMOR_ID
      def idNormal = row.NORMAL_ID
      def target = row.TARGET
      def bamTumor = returnFile(row.TUMOR_BAM)
      // check if using bamTumor.bai or bamTumor.bam.bai
      def baiTumor = returnFile(validateBamIndexFormat(row.TUMOR_BAM))
      // def sizeTumorBamFile = tumorBamFile.size()
      def bamNormal = returnFile(row.NORMAL_BAM)
      def baiNormal = returnFile(validateBamIndexFormat(row.NORMAL_BAM))
      // def sizeNormalBamFile = normalBamFile.size()

      [idTumor, idNormal, target, file(bamTumor), file(baiTumor), file(bamNormal), file(baiNormal)]
    }
  }
  
  
  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
  }

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
  }

  static def check_for_duplicated_rows(pairingFilePath) {
    def entries = []
    file( pairingFilePath ).eachLine { line ->
      if (!line.isEmpty()){
        entries << line
      }
    }
    return entries.toSet().size() == entries.size()
  }


  // check fileIDs are unique in input mapping *tsv
  // not functional for now
  static def checkForUniqueSampleLanes(inputFilename) {
    def totalList = []
    // parse tsv
    file(inputFilename).eachLine { line ->
        if (!line.isEmpty()) {
            def (sample, target, fastqpe1, fastqpe2) = line.split(/\t/)
            totalList << sample + "_" + file(fastqpe1).baseName
        }
    }
    // remove header 'SAMPLE_LANE'
    totalList.removeAll{ it == 'SAMPLE_LANE'} 
    return totalList.size() == totalList.unique().size()
  }

}
