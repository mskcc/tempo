#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - dellyCall
 - dellyFilter
 - CreateIntervalBeds
 - runMutect2
 - indexVCF
 - runMutect2Filter
 - combineMutect2VCF
 - runManta
 - runStrelka
 - combineChannel
 - runBCFToolsFilterNorm
 - runBCFToolsMerge
 - runVCF2MAF
 - doSNPPileup
 - doFacets
 - runMsiSensor
*/


/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

referenceMap = defineReferenceMap()

bamFiles = Channel.empty()

tsvFile = file(tsvPath)

bamFiles = extractBamFiles(tsvFile)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

// --- Run Mutect2
(sampleIdsForIntervalBeds, bamFiles) = bamFiles.into(2)

process CreateScatteredIntervals {

  publishDir "${params.outDir}/intervals", mode: params.publishDirMode

  input:
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
      ])
    set file(idtTargets), file(agilentTargets), file(wgsIntervals) from Channel.value([
      referenceMap.idtTargets,
      referenceMap.agilentTargets,
      referenceMap.wgsTargets
      ])
    set file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervalsIndex) from Channel.value([
      referenceMap.idtTargetsIndex,
      referenceMap.agilentTargetsIndex,
      referenceMap.wgsTargetsIndex
      ])

  output:
    file("agilent*.interval_list") into agilentIntervals mode flatten
    file("idt*.interval_list") into idtIntervals mode flatten
    file("wgs*.interval_list") into wgsIntervals mode flatten

  when: "mutect2" in tools

  script:
  """
  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${agilentTargets} \
    --scatter-count 3 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output agilent

  for i in agilent/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i agilent-\$BASENAME
  done

  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${idtTargets} \
    --scatter-count 3 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output idt

  for i in idt/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i idt-\$BASENAME
  done

  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${wgsIntervals} \
    --scatter-count 3 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output wgs 

  for i in wgs/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i wgs-\$BASENAME
  done
  """
}

( bamsForMutect2Intervals, bamFiles ) = bamFiles.into(2)

//Flipping assay and target in bamList so that it can be combined with the intervals
//Could probably get away with leaving it second and combining, below, by: 1
bamList = bamsForMutect2Intervals
  .map { row ->
    def assay = row[0]
    def target = row[1]
    def idTumor = row[2]
    def idNormal = row[3]
    def bamTumor = returnFile(row[4])
    def bamNormal = returnFile(row[5])
    def baiTumor = returnFile(row[6])
    def baiNormal = returnFile(row[7])
    [ target, assay, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
agilentIList = agilentIntervals.map{ n -> [ "agilent", n ] }
idtIList = idtIntervals.map{ n -> [ "idt", n ] }
wgsIList = wgsIntervals.map{ n -> [ "wgs", n ] }

( aBamList, iBamList, wBamList ) = bamList.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 0).unique() 
bMergedChannel = iBamList.combine(idtIList, by: 0).unique() 
wMergedChannel = wBamList.combine(wgsIList, by: 0).unique() 

mergedChannel = aMergedChannel.concat( bMergedChannel, wMergedChannel)

process RunMutect2 {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2", mode: params.publishDirMode

  input:
    set target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannel 
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz"), file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz.tbi") into mutect2Output

  when: 'mutect2' in tools

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  # Wrong intervals set here
  gatk --java-options -Xmx8g \
    Mutect2 \
    --reference ${genomeFile} \
    --intervals ${intervalBed} \
    --input ${bamTumor} \
    --tumor-sample ${idTumor} \
    --input ${bamNormal} \
    --normal-sample ${idNormal} \
    --output ${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz
  """
}

process RunMutect2Filter {
  tag {idTumor + "_vs_" + idNormal + '_' + mutect2Vcf.baseName}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2_filter", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex) from mutect2Output

  output:
    set idTumor, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forMutect2Combine mode flatten

  when: 'mutect2' in tools && "filter" in tools

  script:
  prefix = "${mutect2Vcf}".replaceFirst('.vcf.gz', '')
  """
  gatk --java-options -Xmx8g \
    FilterMutectCalls \
    --variant ${mutect2Vcf} \
    --stats ${prefix}.Mutect2FilteringStats.tsv \
    --output ${prefix}.filtered.vcf.gz
  """
}

forMutect2Combine = forMutect2Combine.groupTuple(by: [0,1,2])

process CombineMutect2Vcf {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcfOutput

  when: 'mutect2' in tools && 'combine' in tools

  outfile="${idTumor}_vs_${idNormal}.mutect2.filtered.vcf.gz"

  script:
  """
  bcftools concat \
    --allow-overlaps \
    ${mutect2Vcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all | \
  bcftools view \
    --samples ${idNormal},${idTumor} \
    --output-type z \
    --output-file ${outfile}

  tabix --preset vcf ${outfile}
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
  result_array = [
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
    'msiSensorList'    : checkParamReturnFile("msiSensorList"),
    'svCallingExcludeRegions' : checkParamReturnFile("svCallingExcludeRegions"),
    'svCallingIncludeRegions' : checkParamReturnFile("svCallingIncludeRegions"),
    'svCallingIncludeRegionsIndex' : checkParamReturnFile("svCallingIncludeRegionsIndex"),
    'idtTargets' : checkParamReturnFile("idtTargets"),
    'idtTargetsIndex' : checkParamReturnFile("idtTargetsIndex"),
    'agilentTargets' : checkParamReturnFile("agilentTargets"),
    'agilentTargetsIndex' : checkParamReturnFile("agilentTargetsIndex"),
    'wgsTargets' : checkParamReturnFile("wgsTargets"),
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex")
  ]

  if (!params.test) {
    result_array << ['vepCache'                 : checkParamReturnFile("vepCache")]
    // for SNP Pileup
    result_array << ['facetsVcf'        : checkParamReturnFile("facetsVcf")]
    // intervals file for spread-and-gather processes
    result_array << ['intervals'        : checkParamReturnFile("intervals")]
    // files for CombineChannel, needed by bcftools annotate
    result_array << ['repeatMasker'    : checkParamReturnFile("repeatMasker")]
    result_array << ['repeatMaskerIndex'    : checkParamReturnFile("repeatMaskerIndex")]
    result_array << ['mapabilityBlacklist' : checkParamReturnFile("mapabilityBlacklist")]
    result_array << ['mapabilityBlacklistIndex' : checkParamReturnFile("mapabilityBlacklistIndex")]
    // isoforms needed by vcf2maf
    result_array << ['isoforms' : checkParamReturnFile("isoforms")]
    // PON files
    result_array << ['exomePoN' : checkParamReturnFile("exomePoN")]
    result_array << ['exomePoNIndex' : checkParamReturnFile("exomePoNIndex")]
    result_array << ['wgsPoN' : checkParamReturnFile("wgsPoN")]
    result_array << ['wgsPoNIndex' : checkParamReturnFile("wgsPoNIndex")]
  }
  return result_array
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "assay targets idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    checkNumberOfItem(row, 8)
    def assay = row.ASSAY
    def target = row.TARGET
    def idTumor = row.TUMOR_ID
    def idNormal = row.NORMAL_ID
    def bamTumor = returnFile(row.TUMOR_BAM)
    def bamNormal = returnFile(row.NORMAL_BAM)
    def baiTumor = returnFile(row.TUMOR_BAI)
    def baiNormal = returnFile(row.NORMAL_BAI)
    checkFileExtension(bamTumor,".bam")
    checkFileExtension(bamNormal,".bam")
    checkFileExtension(baiTumor,".bai")
    checkFileExtension(baiNormal,".bai")
    [ assay, target, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
}

// Check file extension
def checkFileExtension(it, extension) {
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
  if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Return file if it exists
def returnFile(it) {
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}
