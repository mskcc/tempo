#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - doSNPPileup
 - doFacets
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

process doSNPPileup {
  tag { "SNPPileup_" + idTumor + "_" + idNormal }   // The tag directive allows you to associate each process executions with a custom label

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamFiles
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set idTumor, idNormal, file("${output_filename}") into SNPPileup

  script:
  output_filename = idTumor + "_" + idNormal + ".snppileup.dat.gz"
  """
  snp-pileup -A -P 50 --gzip "${facetsVcf}" "${output_filename}" "${bamTumor}" "${bamNormal}"
  """
}

process doFacets {

  publishDir "${ params.outDir }"

  input:
    set idTumor, idNormal, file("${idTumor}_${idNormal}.snppileup.dat.gz") from SNPPileup

  output:
    file("*.*") into FacetsOutput

  script:
  snp_pileup_prefix = idTumor + "_" + idNormal
  counts_file = "${snp_pileup_prefix}.snppileup.dat.gz"
  genome_value = "hg19"
  TAG = "${snp_pileup_prefix}"
  directory = "."
  """
  /usr/bin/facets-suite/doFacets.R \
  --cval 100 \
  --snp_nbhd 250 \
  --ndepth 35 \
  --min_nhet 25 \
  --purity_cval 500 \
  --purity_snp_nbhd 250 \
  --purity_ndepth 35 \
  --purity_min_nhet 25 \
  --genome "${genome_value}" \
  --counts_file "${counts_file}" \
  --TAG "${TAG}" \
  --directory "${directory}" \
  --R_lib latest \
  --single_chrom F \
  --ggplot2 T \
  --seed 1000 \
  --tumor_id "${idTumor}"
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
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels'      : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
    // for SNP Pileup
    'facetsVcf'        : checkParamReturnFile("facetsVcf") 
  ]
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 6)
    def idTumor = row[0]
    def idNormal = row[1]
    def bamTumor = SarekUtils.returnFile(row[2])
    def bamNormal = SarekUtils.returnFile(row[3])
    def baiTumor = SarekUtils.returnFile(row[4])
    def baiNormal = SarekUtils.returnFile(row[5])

    SarekUtils.checkFileExtension(bamTumor,".bam")
    SarekUtils.checkFileExtension(bamNormal,".bam")
    SarekUtils.checkFileExtension(baiTumor,".bai")
    SarekUtils.checkFileExtension(baiNormal,".bai")

    [ idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
}
