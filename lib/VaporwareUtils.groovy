import static nextflow.Nextflow.file
import nextflow.Channel

class VaporwareUtils {

  static def checkParamReturnFile(item) {
    params."${item}" = params.genomes[params.genome]."${item}"
    return file(params."${item}")
  }

  static def defineReferenceMap() {
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

  // Check which format of BAM index used, input 'it' as BAM file 'bamTumor.bam'
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