#!/usr/bin/env nextflow

/*
=== SOMATIC VARIANT CALLING ===
Processes overview
 - dellyCall
 - dellyFilter
 - makeSamplesFile
 - runMutect2
 - runManta
 - runStrelka
 - doSNPPileup
 - doFacets
 - runMsiSensor 
*/

tsvPath = ''
if (params.sample) tsvPath = params.sample

referenceMap = defineReferenceMap()

bamFiles = Channel.empty()

tsvFile = file(tsvPath)

bamFiles = extractBamFiles(tsvFile)

(bamsForDelly, bamsForMutect2, bamsForManta, bamsForStrelka, bamFilesForSNPPileup, bamsForMakingSampleFile, bamsForMsiSensor, bamsForLumpy) = bamFiles.into(8)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

// --- Run Delly scattered across SV types ---
// --- Runs in two steps, first Call on tumor and normal then filter tumor against normal
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")

process dellyCall {
  tag { "DellyCall_${svType}_" + idTumor + "_vs_" + idNormal }

  publishDir "${params.outDir}/somatic_variants/${idTumor}_vs_${idNormal}/delly"

  input:
    each svType from svTypes
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(dellyExcludeRegions) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.dellyExcludeRegions
    ])

  output:
    set file("${idTumor}_vs_${idNormal}_${svType}.bcf"), file("${idTumor}_vs_${idNormal}_${svType}.bcf.csi"), svType into dellyCallOutput

  when: 'delly' in tools

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${dellyExcludeRegions} \
    --outfile ${idTumor}_vs_${idNormal}_${svType}.bcf \
    ${bamTumor} ${bamNormal}
  """
}

process makeSamplesFile {
  tag { "samplesFile_" + idTumor + "_vs_" + idNormal }

  input: 
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMakingSampleFile 

  output:
    set idTumor, idNormal, file("samples.tsv") into sampleTSVFile

  when: 'delly' in tools

  script:
  """
  echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
  """
} 

process dellyFilter {
  tag { "DellyFilter_${svType}_" + idTumor + "_vs_" + idNormal }

  publishDir "${params.outDir}/somatic_variants/${idTumor}_vs_${idNormal}/delly"

  input:
    set idTumor, idNormal, file("samples.tsv") from sampleTSVFile 
    set file("${idTumor}_vs_${idNormal}_${svType}.bcf"), file("${idTumor}_vs_${idNormal}_${svType}.bcf.csi"), svType from dellyCallOutput

  output:
    set file("${idTumor}_vs_${idNormal}_${svType}.filter.bcf"), file("${idTumor}_vs_${idNormal}_${svType}.filter.bcf.csi") into dellyFilterOutput 

  when: 'delly' in tools

  script:
  """
  delly filter \
    --filter somatic \
    --samples samples.tsv \
    --outfile ${idTumor}_vs_${idNormal}_${svType}.filter.bcf \
    ${idTumor}_vs_${idNormal}_${svType}.bcf
  """
}

// ---------------------- Run MuTect2 

process runMutect2 {
  tag {"MUTECT2_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/mutect2"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMutect2
    set file(genomeFile), file(genomeIndex), file(intervals) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.intervals
    ])

  output:
    set idNormal, idTumor, file("${idTumor}_vs_${idNormal}_somatic.vcf.gz") into mutect2Output

  when: 'mutect2' in tools

  // insert right call regions below

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options -Xmx8g \
    Mutect2 \
    --reference ${genomeFile} \
    --input ${bamTumor} \
    --tumor-sample ${idTumor} \
    --input ${bamNormal} \
    --normal-sample ${idNormal} \
    --output ${idTumor}_vs_${idNormal}_somatic.vcf.gz
  """
}

process indexVCF {
  tag {"INDEXVCF_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/index_vcf"

  input:
    set idTumor, idNormal, file(mutect2Vcf) from mutect2Output

  output:
    set idTumor, idNormal, file(mutect2Vcf), file("${mutect2Vcf.baseName}.gz.tbi") into mutect2IndexedOutput

  when: 'mutect2' in tools

  script:
  """
  tabix -p vcf ${mutect2Vcf} 
  """
}

process runMutect2Filter {
  tag {"MUTECT2FILTER_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/mutect2_filter"

  input:
    set idTumor, idNormal, file(mutect2Vcf), file(mutect2VcfIndex) from mutect2IndexedOutput

  output:
    file("*somatic.filtered.vcf*") into mutect2FilteredOutput

  when: 'mutect2' in tools

  outfile = "${idTumor}_vs_${idNormal}_somatic.filtered.vcf.gz"

  // this process also creates a *.tsv file that you can place write to any path with --stats argument

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    FilterMutectCalls \
    --variant "${mutect2Vcf}" \
    --output "${outfile}" 
  """
}

// ---------------------- Run Manta and Strelka

process runManta {
  tag {"RUNMANTA_" + idTumor + "_" + idNormal}

  publishDir "${params.outDir}/VariantCalling/${idTumor}_${idNormal}/Manta"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
    set file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: 'manta' in tools

  // fix call regions
  // flag with --exome if exome
  script:
  """
  configManta.py \
    --referenceFasta ${genomeFile} \
    --callRegions ${call_regions} \
    --normalBam ${bamNormal} \
    --tumorBam ${bamTumor} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus} \
    --memGb ${task.memory}

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

  publishDir "${params.outDir}/VariantCalling/${idTumor}_${idNormal}/Strelka"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelka
    set file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idNormal, idTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  when: 'manta' in tools && 'strelka2' in tools
  
  // flag with --exome if exome

  script:
  """
  configureStrelkaSomaticWorkflow.py \
    --referenceFasta ${genomeFile} \
    --callRegions ${call_regions} \
    --indelCandidates ${mantaCSI} \
    --tumorBam ${bamTumor} \
    --normalBam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus} \
    --memGb ${task.memory}

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

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/snppileup"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamFilesForSNPPileup
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set sequenceType, idTumor, idNormal, file("${output_filename}") into SNPPileup

  when: 'facets' in tools

  script:
  output_filename = idTumor + "_" + idNormal + ".snppileup.dat.gz"
  """
  snp-pileup \
    --count-orhpans \
    --pseudo-snps 50 \
    --gzip \
    ${facetsVcf} \
    ${output_filename} \
    ${bamTumor} ${bamNormal}
  """
}

process doFacets {
  tag { "DOFACETS_" + idTumor + "_" + idNormal }  

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/facets"

  input:
    set sequenceType, idTumor, idNormal, file(snpPileupFile) from SNPPileup

  output:
    file("*.*") into FacetsOutput

  when: 'facets' in tools

  script:
  snp_pileup_prefix = idTumor + "_" + idNormal
  counts_file = "${snpPileupFile}"
  genome_value = "hg19"
  TAG = "${snp_pileup_prefix}"
  """
  /usr/bin/facets-suite/doFacets.R \
    --cval "${params.facets.cval}" \
    --snp_nbhd "${params.facets.snp_nbhd}" \
    --ndepth "${params.facets.ndepth}" \
    --min_nhet "${params.facets.min_nhet}" \
    --purity_cval "${params.facets.purity_cval}" \
    --purity_snp_nbhd "${params.facets.purity_snp_nbhd}" \
    --purity_ndepth "${params.facets.purity_ndepth}" \
    --purity_min_nhet "${params.facets.purity_min_nhet}" \
    --genome "${params.facets.genome}" \
    --counts_file "${counts_file}" \
    --TAG "${TAG}" \
    --directory "${params.facets.directory}" \
    --R_lib "${params.facets.R_lib}" \
    --single_chrom "${params.facets.single_chrom}" \
    --ggplot2 "${params.facets.ggplot2}" \
    --seed "${params.facets.seed}" \
    --tumor_id ${idTumor}
  """

}

process runMsiSensor {
  tag { "MSISENSOR_" + idTumor + "_" + idNormal }  

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/msisensor"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    file("${output_prefix}*") into msiOutput 

  when: "msisensor" in tools

  script:
  output_prefix = "${idTumor}_${idNormal}"
  """
  msisensor msi \
    -d "${msiSensorList}" \
    -t "${bamTumor}" \
    -n "${bamNormal}" \
    -o "${output_prefix}"
  """
}

process runLumpyExpress {
  tag { "LUMPYEXPRESS_" + idTumor + "_" + idNormal }  

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/lumpyexpress"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForLumpy

  output:
    file("*.vcf") into lumpyExpressOutput

  when: 'lumpyexpress' in tools

  script:
  """
  output_filename=${idTumor}_${idNormal}.lumpyexpress.vcf
  lumpyexpress \
    -B ${bamTumor},${bamNormal} \
    -o "\${output_filename}"
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
    'facetsVcf'        : checkParamReturnFile("facetsVcf"),
    // Microsatellite sites for MSIsensor
    'msiSensorList'    : checkParamReturnFile("msiSensorList"),
    // Genomic regions to exclude for Delly
    'dellyExcludeRegions' : checkParamReturnFile("dellyExcludeRegions")
  ]
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 7)
    def sequenceType = row[0]
    def idTumor = row[1]
    def idNormal = row[2]
    def bamTumor = SarekUtils.returnFile(row[3])
    def bamNormal = SarekUtils.returnFile(row[4])
    def baiTumor = SarekUtils.returnFile(row[5])
    def baiNormal = SarekUtils.returnFile(row[6])

    SarekUtils.checkFileExtension(bamTumor,".bam")
    SarekUtils.checkFileExtension(bamNormal,".bam")
    SarekUtils.checkFileExtension(baiTumor,".bai")
    SarekUtils.checkFileExtension(baiNormal,".bai")

    [ sequenceType, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
}
