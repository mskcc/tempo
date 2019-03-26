#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
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

// ---------------------- Run Delly Call and Filter

sv_variants = Channel.from( "DUP", "BND", "DEL", "INS", "INV" )

( bamsForDelly, bamFiles ) = bamFiles.into(2)

process dellyCall {
  tag { "DELLYCALL_${sv_variant}_" + idTumor + "_" + idNormal }

  publishDir "${params.outDir}/VariantCalling/${idTumor}_${idNormal}/delly_call"

  input:
    each sv_variant from sv_variants
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
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
    set idTumor, idNormal, sv_variant, file("${idTumor}_${idNormal}_${sv_variant}.bcf"), file("${idTumor}_${idNormal}_${sv_variant}.bcf.csi") into dellyCallOutput

  when: 'delly' in tools

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

( bamsForMakingSampleFile, bamFiles ) = bamFiles.into(2)

process makeSamplesFile {
  tag { "SAMPLESFILE_" + idTumor + "_" + idNormal }

  input: 
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMakingSampleFile 

  output:
    file("samples.tsv") into sampleTSVFile

  when: 'delly' in tools

  """
    echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv
  """
} 

dellyCallOutput = dellyCallOutput.spread(sampleTSVFile)

process dellyFilter {
  tag {  idTumor + "_" + idNormal +", " + sv_variant }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/delly_filter"

  input:
    set idTumor, idNormal, sv_variant, file(dellyBcf), file(dellyBcfIndex), file(sampleTsv) from dellyCallOutput

  output:
    set file("*.filter.bcf"), file("*.filter.bcf.csi") into dellyFilterOutput

  when: 'delly' in tools

  outfile="${dellyBcf}".replaceFirst(".bcf",".filter.bcf")

  script:
  """
  delly filter \
    -f somatic \
    -o ${outfile} \
    -s ${sampleTsv} \
    ${dellyBcf}
  """
}

// ---------------------- Run MuTect2 

( sampleIdsForIntervalBeds, bamFiles ) = bamFiles.into(2)

process CreateIntervalBeds {
  tag {intervals.fileName}

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/interval_beds"

  input:
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(intervals) from Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict, referenceMap.intervals])
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForIntervalBeds

  output:
    file 'interval_beds/*.interval_list' into bedIntervals mode flatten

  when: "mutect2" in tools

  script:
  """
  gatk SplitIntervals \
    -R ${genomeFile} \
    -L ${intervals} \
    --scatter-count 10 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    -O interval_beds
  """
}

( bamsForMutect2, bamFiles ) = bamFiles.into(2)
bamsForMutect2Intervals = bamsForMutect2.spread(bedIntervals)

if (params.verbose) bamsForMutect2Intervals = bamsForMutect2Intervals.view {
  "BAMs for Mutect2 with Intervals:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[4].fileName}]"
}

process runMutect2 {
  tag {"MUTECT2_" + intervalBed.baseName + "_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/mutect2"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from bamsForMutect2Intervals
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex
    ])

  output:
    set idTumor, idNormal, file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}_somatic.vcf.gz") into mutect2Output

  when: 'mutect2' in tools

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    Mutect2 \
    -R ${genomeFile}\
    -I ${bamTumor}  -tumor ${idTumor} \
    -I ${bamNormal} -normal ${idNormal} \
    -L ${intervalBed} \
    -O "${idTumor}_vs_${idNormal}_${intervalBed.baseName}_somatic.vcf.gz"
  """
}

process indexVCF {
  tag {"INDEXVCF_" + mutect2Vcf.baseName + "_" + idTumor + "_" + idNormal }

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
  tag {"MUTECT2FILTER_" + mutect2Vcf.baseName + "_" + idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/mutect2_filter"

  input:
    set idTumor, idNormal, file(mutect2Vcf), file(mutect2VcfIndex) from mutect2IndexedOutput

  output:
    file("*somatic.filtered.vcf") into mutect2FilteredOutput
    file("*somatic.filtered.vcf.idx") into mutect2FilteredOutputIndex

  when: 'mutect2' in tools

  outfile="${mutect2Vcf.fileName}".replaceFirst("vcf.gz","filtered.vcf")

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  gatk --java-options "-Xmx8g" \
    FilterMutectCalls \
    --variant "${mutect2Vcf}" \
    --output "${outfile}" 
  """
}

( sampleIdsForMutect2Combine, bamFiles ) = bamFiles.into(2)

process combineMutect2VCF {
  tag {"MUTECT2COMBINE_" + idTumor + "_" + idNormal }

  publishDir "${params.outDir}/VariantCalling/${idTumor}_${idNormal}/mutect2_combined"

  input:
    file(mutect2Vcfs) from mutect2FilteredOutput.collect()
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForMutect2Combine

  output:
    file("${outfile}") into mutect2CombinedVcfOutput

  when: 'mutect2' in tools

  outfile="${idTumor}_${idNormal}.mutect2.filtered.combined.vcf.gz"

  script:
  """
  bcftools concat ${mutect2Vcfs} | bcftools sort --output-type z --output-file ${outfile}
  """
}

// ---------------------- Run Manta and Strelka
( bamsForManta, bamsForStrelka, bamFiles ) = bamFiles.into(3)

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
    set idTumor, idNormal, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
    set file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: 'manta' in tools

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
    set file("*indels.vcf.gz"), file("*indels.vcf.gz.tbi") into strelkaOutputIndels
    set file("*snvs.vcf.gz"), file("*snvs.vcf.gz.tbi") into strelkaOutputSNVs

  when: 'manta' in tools && 'strelka2' in tools

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

// ---------------------- Run bcftools filter, norm, merge

( sampleIdsForCombineChannel, bamFiles ) = bamFiles.into(2)

process combineChannel {
  tag { idTumor + "_" + idNormal }

  input:
    file(mutect2combinedVCF) from mutect2CombinedVcfOutput
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForCombineChannel
    set file(strelkaIndels), file(strelkaIndelsTBI) from strelkaOutputIndels
    set file(strelkaSNV), file(strelkaSNVTBI) from strelkaOutputSNVs

  output:
    set file(mutect2combinedVCF), file(strelkaIndels), file(strelkaSNV) into vcfOutputSet

  when: 'manta' in tools && 'strelka2' in tools && 'mutect2' in tools

  script:
  """
  echo 'placeholder process to make a channel containing vcf data'
  """
}

( sampleIdsForBCFToolsFilterNorm, sampleIdsForBCFToolsMerge, bamFiles ) = bamFiles.into(3)

process runBCFToolsFilterNorm {
  tag { idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/vcf_output"

  input:
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForBCFToolsFilterNorm
    each file(vcf) from vcfOutputSet.flatten()
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    file("*filtered.norm.vcf.gz") into vcfFilterNormOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools

  outfile="${vcf}".replaceFirst('vcf.gz', 'filtered.norm.vcf.gz')

  script:
  """
  tabix -p vcf ${vcf}

  bcftools filter \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --output-type z \
    "${vcf}" | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --output-type z \
    --output "${outfile}" 
  """
}

process runBCFToolsMerge {
  tag { idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/vcf_merged_output"

  input:
    file('*.vcf.gz') from vcfFilterNormOutput.collect()
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForBCFToolsMerge

  output:
    file("*filtered.norm.merge.vcf") into vcfMergedOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools

  script:
  """
  for f in *.vcf.gz
  do
    tabix -p vcf \$f
  done

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type v \
    --output "${idTumor}_${idNormal}.mutect2.strelka2.filtered.norm.merge.vcf" \
    *.vcf.gz
  """
}

( sampleIdsForVCF2MAF, bamFiles ) = bamFiles.into(2)

process runVCF2MAF {
  tag { idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/VariantCalling/${idTumor}_${idNormal}/vcf2maf"

  input:
    file(vcfMerged) from vcfMergedOutput
    set sequenceType, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForVCF2MAF
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vcf2mafFilterVcf), file(vcf2mafFilterVcfIndex), file(vepCache) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.vcf2mafFilterVcf,
      referenceMap.vcf2mafFilterVcfIndex,
      referenceMap.vepCache
    ])

  output:
    file("*.maf") into mafFile

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools

  outfile="${vcfMerged}".replaceFirst(".vcf", ".maf")

  script:
  """
  perl /usr/bin/vcf2maf/vcf2maf.pl \
    --input-vcf ${vcfMerged} \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --filter-vcf ${vcf2mafFilterVcf} \
    --output-maf ${outfile} \
    --ref-fasta ${genomeFile}
  """
}

// ---------------------- Run SNPPileup into doFacets
( bamFilesForSNPPileup, bamFiles ) = bamFiles.into(2)
 
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
  snp-pileup -A -P 50 --gzip "${facetsVcf}" "${output_filename}" "${bamTumor}" "${bamNormal}"
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
  --tumor_id "${idTumor}"
  """
}

( bamsForMsiSensor, bamFiles ) = bamFiles.into(2)

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
  msisensor msi -d "${msiSensorList}" -t "${bamTumor}" -n "${bamNormal}" -o "${output_prefix}"
  """
}

( bamsForLumpy, bamFiles ) = bamFiles.into(2)

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
    // MSI Sensor
    'msiSensorList'    : checkParamReturnFile("msiSensorList"),
    // vcf2maf filter vcf
    'vcf2mafFilterVcf'         : checkParamReturnFile("vcf2mafFilterVcf"),
    'vcf2mafFilterVcfIndex'    : checkParamReturnFile("vcf2mafFilterVcfIndex"),
    'vepCache'                 : checkParamReturnFile("vepCache")
  ]
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "sequenceType idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
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
