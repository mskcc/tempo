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

// --- Run Delly
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDelly, bamFiles) = bamFiles.into(2)

process DellyCall {
  tag {idTumor + "_vs_" + idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/delly"

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, svType, file("${idTumor}_vs_${idNormal}_${svType}.bcf"), file("${idTumor}_vs_${idNormal}_${svType}.bcf.csi") into dellyCallOutput

  when: 'delly' in tools

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idTumor}_vs_${idNormal}_${svType}.bcf \
    ${bamTumor} ${bamNormal}
  """
}

process DellyFilter {
  tag {idTumor + "_vs_" + idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/delly"

  input:
    set idTumor, idNormal, svType, file(dellyBcf), file(dellyBcfIndex) from dellyCallOutput


  output:
    file("*.filter.bcf") into dellyFilterOutput
    file("*.filter.bcf.csi") into dellyFilterIndexedOutput

  when: 'delly' in tools

  outfile="${dellyBcf}".replaceFirst(".bcf",".filter.bcf")

  script:
  """
  echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv

  delly filter \
    --filter somatic \
    --samples samples.tsv \
    --outfile ${outfile} \
    ${dellyBcf}
  """
}

// --- Run Mutect2
(sampleIdsForIntervalBeds, bamFiles) = bamFiles.into(2)

process CreateScatteredIntervals {
  tag {idTumor + "_vs_" + idNormal}

  // publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/intervals"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForIntervalBeds
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
    set file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervals) from Channel.value([
      referenceMap.idtTargetsIndex,
      referenceMap.agilentTargetsIndex,
      referenceMap.wgsTargetsIndex
      ])

  output:
    file('intervals/*.interval_list') into bedIntervals mode flatten

  when: "mutect2" in tools

  script:
  intervals = wgsIntervals
  if(params.exome) {
    if(target == 'agilent') intervals = agilentTargets
    if(target == 'idt') intervals = idtTargets
  }
  """
  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${intervals} \
    --scatter-count 30 \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output intervals
  """
}

(bamsForMutect2, bamFiles) = bamFiles.into(2)
bamsForMutect2Intervals = bamsForMutect2.spread(bedIntervals)

if (params.verbose) bamsForMutect2Intervals = bamsForMutect2Intervals.view {
  "BAMs for Mutect2 with Intervals:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[4].fileName}]"
}

process RunMutect2 {
  tag {idTumor + "_vs_" + idNormal + "_" + intervalBed.baseName}

  // publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from bamsForMutect2Intervals
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz"), file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz.tbi") into mutect2Output

  when: 'mutect2' in tools

  // insert right call regions below

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

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2"

  input:
    set idTumor, idNormal, file(mutect2Vcf), file(mutect2VcfIndex) from mutect2Output

  output:
    file("*filtered.vcf.gz") into mutect2FilteredOutput
    file("*filtered.vcf.gz.tbi") into mutect2FilteredOutputIndex

  when: 'mutect2' in tools

  outfile="${mutect2Vcf}".replaceFirst('vcf.gz', 'filtered.vcf.gz')

  // this process also creates a *.tsv file that you can place write to any path with --stats argument

  script:
  """
  gatk --java-options -Xmx8g \
    FilterMutectCalls \
    --variant ${mutect2Vcf} \
    --output ${outfile}
  """
}

(sampleIdsForMutect2Combine, bamFiles) = bamFiles.into(2)

process CombineMutect2Vcf {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2"

  input:
    file(mutect2Vcf) from mutect2FilteredOutput.collect()
    file(mutect2VcfIndex) from mutect2FilteredOutputIndex.collect()
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForMutect2Combine

  output:
    file("${outfile}") into mutect2CombinedVcfOutput

  when: 'mutect2' in tools

  outfile="${idTumor}_vs_${idNormal}.mutect2.filtered.vcf.gz"

  script:
  """
  # Add norm?
  bcftools concat \
    --allow-overlaps \
    ${mutect2Vcf} | \
    bcftools sort \
    --output-type z \
    --output-file ${outfile}
  """
}

// --- Run Manta
(bamsForManta, bamsForStrelka, bamFiles) = bamFiles.into(3)

process RunManta {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/manta"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions,
      referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    file("*.vcf.gz") into mantaOutput
    file("*.vcf.gz.tbi") into mantaIndexedOutput
    set file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: 'manta' in tools
  script:
  options = ""
  if(params.exome) options = "--exome"

  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --referenceFasta ${genomeFile} \
    --normalBam ${bamNormal} \
    --tumorBam ${bamTumor} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

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

// --- Run Strelka2
process RunStrelka2 {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/strelka2"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelka
    set file(mantaCSI), file(mantaCSIi) from mantaToStrelka
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
    set file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervals) from Channel.value([
      referenceMap.idtTargetsIndex,
      referenceMap.agilentTargetsIndex,
      referenceMap.wgsTargetsIndex
      ])

  output:
    set file("*indels.vcf.gz"), file("*indels.vcf.gz.tbi") into strelkaOutputIndels
    set file("*snvs.vcf.gz"), file("*snvs.vcf.gz.tbi") into strelkaOutputSNVs

  when: 'manta' in tools && 'strelka2' in tools

  script:
  options = "" 
  if(params.exome) options = "--exome"

  intervals = wgsIntervals
  if(params.exome) {
    if(target == 'agilent') intervals = agilentTargets
    if(target == 'idt') intervals = idtTargets
  }
  """
  configureStrelkaSomaticWorkflow.py \
    ${options} \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --indelCandidates ${mantaCSI} \
    --tumorBam ${bamTumor} \
    --normalBam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

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

// --- Process Delly and Manta VCFs 

( sampleIdsForDellyMantaMerge, bamFiles ) = bamFiles.into(2)

process MergeDellyAndManta {
  tag { idTumor + "_" + idNormal }

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/vcf_merged_output"

  input:
    file(dellyFilterData) from dellyFilterOutput.collect()
    file(mantaData) from mantaOutput.collect()
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForDellyMantaMerge

  output:
    file("*filtered.merge.vcf") into vcfDellyMantaMergedOutput

  when: 'manta' in tools && 'delly' in tools

  script:
  """ 
  for f in *.bcf
  do 
    bcftools view --output-type z \$f > \${f%.bcf}.vcf.gz
  done

  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type v \
    --output ${idTumor}_${idNormal}.delly.manta.filtered.merge.vcf \
    *.vcf.gz
  """
}

// --- Process Mutect2 and Strelka2 VCFs
(sampleIdsForCombineChannel, bamFiles) = bamFiles.into(2)

process CombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    file(mutect2combinedVCF) from mutect2CombinedVcfOutput
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForCombineChannel
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

(sampleIdsForBcfToolsFilterNorm, sampleIdsForBcfToolsMerge, bamFiles) = bamFiles.into(3)

process RunBcfToolsFilterNorm {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/vcf_output"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForBcfToolsFilterNorm
    each file(vcf) from vcfOutputSet.flatten()
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    file("*filtered.norm.vcf.gz") into vcfFilterNormOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools

  outfile = "${vcf}".replaceFirst('vcf.gz', 'filtered.norm.vcf.gz')

  script:
  """
  tabix --preset vcf ${vcf}

  bcftools filter \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --output-type z \
    ${vcf} | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --output-type z \
    --output ${outfile}
  """
}

process RunBcfToolsMerge {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/vcf_merged_output"

  input:
    file('*.vcf.gz') from vcfFilterNormOutput.collect()
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForBcfToolsMerge

  output:
    file("*filtered.norm.merge.vcf") into vcfMergedOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools

  script:
  """
  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type v \
    --output ${idTumor}_${idNormal}.mutect2.strelka2.filtered.norm.merge.vcf \
    *.vcf.gz
  """
}

(sampleIdsForVcf2Maf, bamFiles) = bamFiles.into(2)

process RunVcf2Maf {
  tag { idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/${idTumor}_vs_${idNormal}/somatic_variants/mutations"

  input:
    file(vcfMerged) from vcfMergedOutput
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForVcf2Maf
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
  perl /opt/vcf2maf.pl \
    --input-vcf ${vcfMerged} \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data ${vepCache} \
    --filter-vcf ${vcf2mafFilterVcf} \
    --output-maf ${outfile} \
    --ref-fasta ${genomeFile}
  """
}

// --- Run FACETS
(bamFilesForSnpPileup, bamFiles) = bamFiles.into(2)
 
process DoSnpPileup {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamFilesForSnpPileup
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set assay, target, idTumor, idNormal, file("${output_filename}") into SnpPileup

  when: 'facets' in tools

  script:
  output_filename = idTumor + "_" + idNormal + ".snppileup.dat.gz"
  """
  snp-pileup \
    --count-orphans \
    --pseudo-snps 50 \
    --gzip \
    ${facetsVcf} \
    ${output_filename} \
    ${bamTumor} ${bamNormal}
  """
}

process DoFacets {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets"

  input:
    set assay, target, idTumor, idNormal, file(snpPileupFile) from SnpPileup

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

// --- Run MSIsensor
(bamsForMsiSensor, bamFiles) = bamFiles.into(2)

process RunMsiSensor {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/msisensor"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForMsiSensor
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

// --- Run HLA Polysolver 
(bamsForHlaPolysolver, bamFiles) = bamFiles.into(2)

process RunHlaPolysolver {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/hla_polysolver"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForHlaPolysolver

  output:
    file("output/*") into hlaOutput

  when: "hla" in tools
  
  script:
  outDir = "output"
  TMPDIR = "$outDir-nf-scratch"
  """
  cp /home/polysolver/scripts/shell_call_hla_type .
  
  sed -i "171s/TMP_DIR=.*/TMP_DIR=$TMPDIR/" shell_call_hla_type 

  bash shell_call_hla_type \
  ${bamNormal} \
  Unknown \
  1 \
  hg19 \
  STDFQ \
  0 \
  ${outDir} ||  echo "HLA Polysolver did not run successfully and its process has been redirected to generate this file." > ${outDir}/winners.hla.txt 
  """
}

// --- Run Conpair

(bamsForConpair, bamFiles) = bamFiles.into(2)

process RunConpair {
  tag {idTumor + "_vs_" + idNormal}


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
    result_array << ['vcf2mafFilterVcf'         : checkParamReturnFile("vcf2mafFilterVcf")]
    result_array << ['vcf2mafFilterVcfIndex'    : checkParamReturnFile("vcf2mafFilterVcfIndex")]
    result_array << ['vepCache'                 : checkParamReturnFile("vepCache")]
    // for SNP Pileup
    result_array << ['facetsVcf'        : checkParamReturnFile("facetsVcf")]
    // intervals file for spread-and-gather processes
    result_array << ['intervals'        : checkParamReturnFile("intervals")]
  }
  return result_array
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "assay targets idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    SarekUtils.checkNumberOfItem(row, 8)
    def assay = row[0]
    def target = row[1]
    def idTumor = row[2]
    def idNormal = row[3]
    def bamTumor = SarekUtils.returnFile(row[4])
    def bamNormal = SarekUtils.returnFile(row[5])
    def baiTumor = SarekUtils.returnFile(row[6])
    def baiNormal = SarekUtils.returnFile(row[7])

    SarekUtils.checkFileExtension(bamTumor,".bam")
    SarekUtils.checkFileExtension(bamNormal,".bam")
    SarekUtils.checkFileExtension(baiTumor,".bai")
    SarekUtils.checkFileExtension(baiNormal,".bai")

    [ assay, target, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
  }
}
