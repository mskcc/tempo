#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - dellyCall
 - dellyFilter
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

( bamsForDelly, bamsForMutect2, bamsForManta, bamsForStrelka, bamFilesForSNPPileup, bamsForMakingSampleFile ) = bamFiles.into(6)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// ---------------------- Run Delly Call and Filter

sv_variants = Channel.from( "DUP", "BND", "DEL", "INS", "INV" )

process dellyCall {
  tag { "DELLYCALL_${sv_variant}_" + idTumor + "_" + idNormal }

  publishDir "${params.outDir}/VariantCalling/delly_call"

  input:
    each sv_variant from sv_variants
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.knownIndels,
      referenceMap.knownIndelsIndex,
      referenceMap.intervals,
    ])

  output:
    set file("${idTumor}_${idNormal}_${sv_variant}.bcf"), file("${idTumor}_${idNormal}_${sv_variant}.bcf.csi"), sv_variant into dellyCallOutput

  """
  outfile="${idTumor}_${idNormal}_${sv_variant}.bcf" 
  delly call \
    -t "${sv_variant}" \
    -o "\${outfile}" \
    -g ${genomeFile} \
    ${bamTumor} \
    ${bamNormal}
  """
}

process makeSamplesFile {
  tag { "SAMPLESFILE_" + idTumor + "_" + idNormal }

  input: 
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMakingSampleFile 

  output:
    set idTumor, idNormal, file("samples.tsv") into sampleTSVFile

  """
    echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
  """
} 

process dellyFilter {
  tag { "DELLYFILTER_${sv_variant}_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/delly_filter"

  input:
    set idTumor, idNormal, file("samples.tsv") from sampleTSVFile 
    set file("${idTumor}_${idNormal}_${sv_variant}.bcf"), file("${idTumor}_${idNormal}_${sv_variant}.bcf.csi"), sv_variant from dellyCallOutput

  output:
    set file("${idTumor}_${idNormal}_${sv_variant}.filter.bcf"), file("${idTumor}_${idNormal}_${sv_variant}.filter.bcf.csi") into dellyFilterOutput 

  """
  delly_call_file="${idTumor}_${idNormal}_${sv_variant}.bcf" 
  outfile="${idTumor}_${idNormal}_${sv_variant}.filter.bcf" 
  delly filter \
    -f somatic \
    -o "\${outfile}" \
    -s "samples.tsv" \
    "\${delly_call_file}"
  """
}

// ---------------------- Run MuTect2 

process runMutect2 {
  tag {"MUTECT2_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/mutect2"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMutect2
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set idNormal, idTumor, file("${idTumor}_vs_${idNormal}_somatic.vcf") into mutect2Output

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    Mutect2 \
    -R ${genomeFile}\
    -I ${bamTumor}  -tumor ${idTumor} \
    -I ${bamNormal} -normal ${idNormal} \
    -O "${idTumor}_vs_${idNormal}_somatic.vcf"
  """
}

// ---------------------- Run Manta and Strelka

process runManta {
  tag {"RUNMANTA_" + idTumor + "_" + idNormal}

  publishDir "${params.outDir}/VariantCalling/Manta"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
    set file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  script:
  """
  configManta.py \
  --normalBam ${bamNormal} \
  --tumorBam ${bamTumor} \
  --reference ${genomeFile} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    Manta_${idTumor}_vs_${idNormal}.somaticSV.vcf.gz.tbi
  """
}

process runStrelka {
  tag {"RUNSTRELKA_" + idTumor + "_" + idNormal}

  publishDir "${params.outDir}/VariantCalling/Strelka"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelka
    set file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  script:
  """
  configureStrelkaSomaticWorkflow.py \
  --tumor ${bamTumor} \
  --normal ${bamNormal} \
  --referenceFasta ${genomeFile} \
  --indelCandidates ${mantaCSI} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz.tbi
  """
}

// ---------------------- Run SNPPileup into doFacets

process doSNPPileup {
  tag { "SNPPILEUP_" + idTumor + "_" + idNormal }  

  publishDir "${ params.outDir }/VariantCalling/snppileup"

  input:
    set idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamFilesForSNPPileup
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
  tag { "DOFACETS_" + idTumor + "_" + idNormal }  

  publishDir "${ params.outDir }/VariantCalling/facets"

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
