if( awscli ) {
    import static nextflow.Nextflow.file
import nextflow.Channel

class VaporwareUtils {

  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
  }

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
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

  // Compare each parameter with a list of parameters
  static def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Return element in list of allowed params
  static def checkParams(it) {
    return it in [
      'annotate-tools',
      'annotate-VCF',
      'annotateTools',
      'annotateVCF',
      'awsqueue',
      'awsqueue_tiny',
      'build',
      'call-name',
      'callName',
      'contact-mail',
      'contactMail',
      'container-path',
      'containerPath',
      'containers',
      'docker',
      'download',
      'explicit-bqsr-needed',
      'explicitBqsrNeeded',
      'genome_base',
      'genome',
      'genomes',
      'help',
      'localReportDir',
      'local-report-dir',
      'markdup_java_options',
      'max_cpus',
      'max_memory',
      'max_time',
      'more',
      'nf-required-version',
      'nfRequiredVersion',
      'no-BAMQC',
      'no-GVCF',
      'no-reports',
      'noBAMQC',
      'noGVCF',
      'noReports',
      'nucleotides-per-second',
      'nucleotidesPerSecond',
      'only-QC',
      'onlyQC',
      'out-dir',
      'outDir',
      'params',
      'project',
      'publish-dir-mode',
      'publishDirMode',
      'push',
      'ref-dir',
      'refDir',
      'repository',
      'run-time',
      'runTime',
      'sample-dir',
      'sample',
      'sampleDir',
      'sequencing_center',
      'single-CPUMem',
      'singleCPUMem',
      'singularity',
      'step',
      'strelka-BP',
      'strelkaBP',
      'tag',
      'target-BED',
      'targetBED',
      'test',
      'tools',
      'total-memory',
      'totalMemory',
      'vcflist',
      'verbose',
      'version']
  }

  // Loop through all the references files to check their existence
  static def checkReferenceMap(referenceMap) {
    referenceMap.every {
      referenceFile, fileToCheck ->
      VaporwareUtils.checkRefExistence(referenceFile, fileToCheck)
    }
  }

  // Loop through all the references files to check their existence

  static def checkRefExistence(referenceFile, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ VaporwareUtils.checkRefExistence(referenceFile, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
			println  "Missing references: ${referenceFile} ${fileToCheck}"
      return false
    }
    return true
  }

  // Define map of directories
  // Needs to be changed
  static def defineDirectoryMap(outDir) {
    return [
    'duplicateMarked'  : "${outDir}/Preprocessing/DuplicateMarked",
    'recalibrated'     : "${outDir}/Preprocessing/Recalibrated",
    'ascat'            : "${outDir}/VariantCalling/Ascat",
    'freebayes'        : "${outDir}/VariantCalling/FreeBayes",
    'gvcf-hc'          : "${outDir}/VariantCalling/HaplotypeCallerGVCF",
    'haplotypecaller'  : "${outDir}/VariantCalling/HaplotypeCaller",
    'manta'            : "${outDir}/VariantCalling/Manta",
    'mutect2'          : "${outDir}/VariantCalling/MuTect2",
    'strelka'          : "${outDir}/VariantCalling/Strelka",
    'strelkabp'        : "${outDir}/VariantCalling/StrelkaBP",
    'snpeff'           : "${outDir}/Annotation/SnpEff",
    'vep'              : "${outDir}/Annotation/VEP",
    'bamQC'            : "${outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${outDir}/Reports/BCFToolsStats",
    'fastQC'           : "${outDir}/Reports/FastQC",
    'markDuplicatesQC' : "${outDir}/Reports/MarkDuplicates",
    'multiQC'          : "${outDir}/Reports/MultiQC",
    'samtoolsStats'    : "${outDir}/Reports/SamToolsStats",
    'snpeffReports'    : "${outDir}/Reports/SnpEff",
    'vcftools'         : "${outDir}/Reports/VCFTools",
    'version'          : "${outDir}/Reports/ToolsVersion"
    ]
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  static def extractBams(tsvFile, mode) {
    Channel.from(tsvFile)
      .splitCsv(sep: '\t')
      .map { row ->
        VaporwareUtils.checkNumberOfItem(row, 6)
        def idPatient = row[0]
        def gender    = row[1]
        def status    = VaporwareUtils.returnStatus(row[2].toInteger())
        def idSample  = row[3]
        def bamFile   = VaporwareUtils.returnFile(row[4])
        def baiFile   = VaporwareUtils.returnFile(row[5])

        VaporwareUtils.checkFileExtension(bamFile,".bam")
        VaporwareUtils.checkFileExtension(baiFile,".bai")

        if (mode == "germline") return [ idPatient, status, idSample, bamFile, baiFile ]
        else return [ idPatient, gender, status, idSample, bamFile, baiFile ]
      }
  }

  // Compare params to list of verified params
  static def isAllowedParams(params) {
    def test = true
    params.each{
      if (!checkParams(it.toString().split('=')[0])) {
        println "params ${it.toString().split('=')[0]} is unknown"
        test = false
      }
    }
    return test
  }

  // Return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
  }

  // Return status [0,1]
  // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }

}
