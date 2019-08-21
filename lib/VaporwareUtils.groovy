import static nextflow.Nextflow.file
import nextflow.Channel

// NOTE: 
// static methods in Groovy are meant to be accessed directly from the class
// e.g. VaporewareUtils.debug()
//

class VaporwareUtils {

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

  static def extractFastq(tsvFile) {
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

  // Check which format of BAM index used, input 'it' as BAM file 'bamTumor.bam'
  // not a static method, as currently written
  static def validateBamIndexFormat(it) {
    bamFilename = it.take(it.lastIndexOf('.'))
    // Check BAM index extension
    if (file(bamFilename + ".bai").exists()){
      return(file("${bamFilename}.bai"))
    } else if (file(bamFilename + ".bam.bai").exists()){
      return(file("${bamFilename}.bam.bai"))
    } else {
      println "ERROR: Cannot find BAM indices for ${it}. Please index BAMs in the same directory with 'samtools index' and re-run the pipeline."
      exit 1
    }
  }

  
  static def extractBAM(tsvFile) {
    Channel.from(tsvFile)
    .splitCsv(sep: '\t', header: true)
    .map { row ->
      checkNumberOfItem(row, 6)
      def idTumor = row.TUMOR_ID
      def idNormal = row.NORMAL_ID
      def assay = row.ASSAY
      def target = row.TARGET
      def bamTumor = returnFile(row.TUMOR_BAM)
      // check if using bamTumor.bai or bamTumor.bam.bai
      def baiTumor = returnFile(validateBamIndexFormat(row.TUMOR_BAM))
      // def sizeTumorBamFile = tumorBamFile.size()
      def bamNormal = returnFile(row.NORMAL_BAM)
      def baiNormal = returnFile(validateBamIndexFormat(row.NORMAL_BAM))
      // def sizeNormalBamFile = normalBamFile.size()

      [assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)]
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

  static def check_for_mixed_assay(mappingFilePath) {
    def wgs = false
    def wes = false
    file( mappingFilePath ).eachLine { line -> 
      if (line.toLowerCase().contains('\tgenome\t') || line.toLowerCase().contains('\twgs\t')) {
        wgs = true
      }
      if (line.toLowerCase().contains('\texome\t') || line.toLowerCase().contains('\twes\t')) {
        wes = true
      }
    return !(wgs && wes)
    }
  }

  // check lane names are unique in input mapping *tsv 
  static def checkForUniqueSampleLanes(inputFilename) {
    def totalList = []
    // parse tsv
    file(inputFilename).eachLine { line ->
        if (!line.isEmpty()){
            def (sample, lane, assay, target, fastqpe1, fastqpe2) = line.split(/\t/)
            totalList << sample + "_" + lane
        }
    }
    // remove header 'SAMPLE_LANE'
    totalList.removeAll{ it == 'SAMPLE_LANE'} 
    return totalList.size() == totalList.unique().size()
  }

}




