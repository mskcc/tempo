
def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  if(params."${item}" == null){println "${item} is not found in reference map"; exit 1}
  if(file(params."${item}", checkIfExists: false) == []){println "${item} is not found; glob pattern produces empty list"; exit 1}
  return file(params."${item}", checkIfExists: true)
}

def defineReferenceMap() {
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
  ]

  if (workflow.profile != "test") {
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
  }
  return result_array
}

def loadTargetReferences(){
  def result_array = [:]
  new File(params.targets_base).eachDir{ i -> 
    def target_id = i.getBaseName()
    result_array["${target_id}"] = [:]
    for ( j in params.targets.keySet()) { // baitsInterval, targetsInterval, targetsBedGz, targetsBedGzTbi, codingBed
      result_array."${target_id}" << [ ("$j".toString()) : evalTargetPath(j,target_id)]
    }
  }
  return result_array
}

def evalTargetPath(item,target_id){
  def templateString = params.targets."${item}"
  if(templateString == null){println "${item} is not found in targets' map"; exit 1}
  def res = evaluate("def targets_id=\"$target_id\" ; template=\"$templateString\"")
  if(file(res, checkIfExists: false) == []){println "${item} is not found; glob pattern produces empty list"; exit 1}
  return file(file(res, checkIfExists: true).toAbsolutePath().toRealPath())
}

def touchInputs() {
  new Timer().schedule({
  def timeNow = new Date().getTime()
  limitInputLines = chunkSizeLimit + ( ((timeNow - startEpoch)/60000) * (chunkSizeLimit / params.touchInputsInterval) )
  for ( i in epochMap.keySet() ){
    fileEpoch = file(i).lastModified()
    if (( fileEpoch > epochMap[i]) || (chunkSizeLimit > 0 )) {
      epochMap[i] = fileEpoch
      "touch -ca ${i}".execute()
    }
  }
} as TimerTask, 1000, params.touchInputsInterval * 60 * 1000 ) // convert minutes to milliseconds
}

def watchMapping(tsvFile, assayType) {
  def index = 1 
  Channel.watchPath( tsvFile, 'create, modify' )
	 .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
	 .view()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: false)
              def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: false)
              def numOfPairs = row.NUM_OF_PAIRS.toInteger()
              if(!TempoUtils.checkTarget(target, assayType, targetsMap.keySet())){}
              if(!TempoUtils.checkNumberOfItem(row, 5, tsvFile)){}

              [idSample, numOfPairs, target, fastqFile1, fastqFile2]
      }
      .map{ idSample, numOfPairs, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, numOfPairs), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchBamMapping(tsvFile, assayType){
  def index = 1
  Channel.watchPath( tsvFile, 'create, modify' )
	 .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def bam = file(row.BAM, checkIfExists: false)
              def bai = file(row.BAI, checkIfExists: false)
              if(!TempoUtils.checkTarget(target, assayType, targetsMap.keySet())){}
              if(!TempoUtils.checkNumberOfItem(row, 4, tsvFile)){}

              [idSample, target, bam, bai]
      }
      .map{ idSample, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, 1), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchPairing(tsvFile){
  Channel.watchPath( tsvFile, 'create, modify' )
	 .splitCsv(sep: '\t', header: true)
	 .unique()
         .map { row ->
              def TUMOR_ID = row.TUMOR_ID
              def NORMAL_ID = row.NORMAL_ID
              if(!TempoUtils.checkNumberOfItem(row, 2, tsvFile)){}

              [TUMOR_ID, NORMAL_ID]
         }
	 .unique()
}

def watchAggregateWithPath(tsvFile) {
  def index = 1 
  Channel.watchPath(tsvFile, 'create, modify')
     .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              def path = row.PATH
              if(!TempoUtils.checkNumberOfItem(row, 5, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal, path]
         }
         .map { cohort, cohortSize, idTumor, idNormal, path
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal, path)
         }
         .transpose()
	 .unique()
}

def watchAggregate(tsvFile) {
  Channel.watchPath(file(runAggregate), 'create, modify')
     .splitCsv(sep: '\t', header: true)
	 .unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              if(!TempoUtils.checkNumberOfItem(row, 4, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal]
         }
         .map { cohort, cohortSize, idTumor, idNormal
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal)
         }
         .transpose()
	 .unique()
}
