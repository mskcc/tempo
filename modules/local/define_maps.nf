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
    result_array << ['annotSVref' : checkParamReturnFile('annotSVref')]
    // splice sites for locating cDNA contam
    result_array << ['spliceSites' : checkParamReturnFile('spliceSites')]
    
  }
  return result_array
}

def loadTargetReferences(){
  def result_array = [:]
  new File(params.targets_base).eachDir{ i -> 
    def target_id = i.getBaseName()
    if (params.assayType == "genome" && target_id != "wgs" ){ return }
    if (params.assayType != "genome" && target_id == "wgs" ){ return }
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
