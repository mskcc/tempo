#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - SomaticDellyCall
 - CreateIntervalBeds
 - RunMutect2
 - RunMutect2Filter
 - SomaticCombineMutect2VCF
 - SomaticRunManta
 - SomaticRunStrelka
 - SomaticCombineChannel
 - SomaticRunBCFToolsFilterNorm
 - SomaticRunBCFToolsMerge
 - SomaticAnnotateMaf
 - SomaticDoSNPPileup
 - DoFacets 
 - RunMsiSensor
 - Polysolver
 - LOHHLA
 - RunConpair
 - RunMutationSignatures
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

// parse --tools parameter for downstream 'when' conditionals, e.g. when: `` 'delly ' in tools
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

if('strelka2' in tools) {
  tools.add('manta')
}

// --- Run Delly

svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDelly, bamFiles) = bamFiles.into(2)

process SomaticDellyCall {
  tag {idTumor + "_vs_" + idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/delly", mode: params.publishDirMode

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}_${svType}.filter.bcf")  into dellyFilterOutput

  when: 'delly' in tools && runSomatic

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idTumor}_vs_${idNormal}_${svType}.bcf \
    ${bamTumor} ${bamNormal}

  echo "${idTumor}\ttumor\n${idNormal}\tcontrol" > samples.tsv

  delly filter \
    --filter somatic \
    --samples samples.tsv \
    --outfile ${idTumor}_vs_${idNormal}_${svType}.filter.bcf \
    ${idTumor}_vs_${idNormal}_${svType}.bcf
  """
}




process CreateScatteredIntervals {

  //publishDir "${params.outDir}/intervals", mode: params.publishDirMode

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
  scatterCount = 10
  """
  gatk SplitIntervals \
    --reference ${genomeFile} \
    --intervals ${agilentTargets} \
    --scatter-count ${scatterCount} \
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
    --scatter-count ${scatterCount} \
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
    --scatter-count ${scatterCount} \
    --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
    --output wgs 

  for i in wgs/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i wgs-\$BASENAME
  done
  """
}


// --- Run Mutect2


(bamsForMutect2Intervals, bamFiles) = bamFiles.into(2)

//Associating interval_list files with BAM files, putting them into one channel
agilentIList = agilentIntervals.map{ n -> [ n, "agilent" ] }
idtIList = idtIntervals.map{ n -> [ n, "idt" ] }
wgsIList = wgsIntervals.map{ n -> [ n, "wgs" ] }

(aBamList, iBamList, wBamList) = bamsForMutect2Intervals.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 1).unique() 
bMergedChannel = iBamList.combine(idtIList, by: 1).unique() 
wMergedChannel = wBamList.combine(wgsIList, by: 1).unique() 

mergedChannel = aMergedChannel.concat( bMergedChannel, wMergedChannel)

process RunMutect2 {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2", mode: params.publishDirMode

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannel 
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forMutect2Combine mode flatten


  when: 'mutect2' in tools

  script:
  mutect2Vcf = "${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz"
  prefix = "${mutect2Vcf}".replaceFirst('.vcf.gz', '')
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
    --output ${mutect2Vcf}

  gatk --java-options -Xmx8g \
    FilterMutectCalls \
    --variant ${mutect2Vcf} \
    --stats ${prefix}.Mutect2FilteringStats.tsv \
    --output ${prefix}.filtered.vcf.gz
  """
}

//Formatting the channel to be keyed by idTumor, idNormal, and target
forMutect2Combine = forMutect2Combine.groupTuple(by: [0,1,2])

process SomaticCombineMutect2Vcf {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcfOutput

  when: 'mutect2' in tools

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

// --- Run Manta
(bamsForManta, bamsForStrelka, bamFiles) = bamFiles.into(3)

process SomaticRunManta {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/manta", mode: params.publishDirMode

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
    set idTumor, idNormal, target, file("*.vcf.gz") into mantaOutput mode flatten
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka mode flatten

  when: 'manta' in tools
  script:
  options = ""
  if(params.assayType == "exome") options = "--exome"

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

mantaOutput = mantaOutput.groupTuple(by: [0,1,2])

dellyFilterOutput = dellyFilterOutput.groupTuple(by: [0,1,2])

dellyMantaCombineChannel = dellyFilterOutput.combine(mantaOutput, by: [0,1,2]).unique()

// --- Process Delly and Manta VCFs 

(sampleIdsForDellyMantaMerge, bamFiles) = bamFiles.into(2)

process SomaticMergeDellyAndManta {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/vcf_merged_output", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(dellyBcfs), file(mantaFile) from dellyMantaCombineChannel

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

// --- Run Strelka2

mantaToStrelka = mantaToStrelka.groupTuple(by: [0,1,2])

process SomaticRunStrelka2 {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/strelka2", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(mantaCSI), file(mantaCSIi) from mantaToStrelka
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
    set idTumor, idNormal, target, file("*indels.vcf.gz"), file("*indels.vcf.gz.tbi"), file("*snvs.vcf.gz"), file("*snvs.vcf.gz.tbi") into strelkaOutput mode flatten

  when: 'manta' in tools && 'strelka2' in tools

  script:
  options = ""
  intervals = wgsIntervals
  if(params.assayType == "exome") {
    options = "--exome"
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

strelkaOutput = strelkaOutput.groupTuple(by: [0,1,2])

// --- Process Mutect2 and Strelka2 VCFs

process SomaticMergeStrelka2Vcfs {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/strelka2", mode: params.publishDirMode

  input: 
    set idTumor, idNormal, target, file(strelkaIndels), file(strelkaIndelsIndex), file(strelkaSNVs), file(strelkaSNVsIndex) from strelkaOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file('*.vcf.gz'), file('*.vcf.gz.tbi') into strelkaOutputMerged

  when: 'manta' in tools && 'strelka2' in tools

  script:
  prefix = "${idTumor}_${idNormal}_${target}.strelka.merged"
  outfile = "${prefix}.filtered.vcf.gz"
  """
  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    ${strelkaIndels} ${strelkaSNVs} | \
  bcftools reheader \
    --samples samples.txt | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --output-type z \
    --output ${outfile}

  tabix --preset vcf ${outfile}
  """
}

mutectStrelkaChannel = mutect2CombinedVcfOutput.combine(strelkaOutputMerged, by: [0,1,2]).unique()

process SomaticCombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutectCombinedVcf), file(mutectCombinedVcfIndex), file(strelkaVcf), file(strelkaVcfIndex) from mutectStrelkaChannel
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker,
      referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist,
      referenceMap.mapabilityBlacklistIndex
    ])
    set file(exomePoN), file(wgsPoN), file(exomePoNIndex), file(wgsPoNIndex) from Channel.value([
      referenceMap.exomePoN,
      referenceMap.wgsPoN,
      referenceMap.exomePoNIndex,
      referenceMap.wgsPoNIndex
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf,
      referenceMap.gnomadWesVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.pass.vcf") into vcfMergedOutput

  when: 'strelka2' in tools && 'mutect2' in tools

  script:
  isec_dir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWesVcf
  if (target != 'genome') {
    pon = exomePoN
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > vcf.rm.header
  echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > vcf.map.header
  echo -e "##INFO=<ID=PoN,Number=1,Type=Integer,Description=\"Count in panel of normals\">" > vcf.pon.header

  bcftools isec \
    --output-type z \
    --prefix ${isec_dir} \
    ${mutectCombinedVcf} ${strelkaVcf}

  bcftools annotate \
    --annotations ${isec_dir}/0003.vcf.gz \
    --include 'FILTER!=\"PASS\"' \
    --mark-sites \"+Strelka2FILTER\" \
    -k \
    --output-type z \
    --output ${isec_dir}/0003.annot.vcf.gz \
    ${isec_dir}/0003.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0000.vcf.gz \
    --mark-sites +MuTect2 \
    --output-type z \
    --output ${isec_dir}/0000.annot.vcf.gz \
    ${isec_dir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0002.vcf.gz \
    --mark-sites \"+MuTect2;Strelka2\" \
    --output-type z \
    --output ${isec_dir}/0002.tmp.vcf.gz \
    ${isec_dir}/0002.vcf.gz

  tabix --preset vcf ${isec_dir}/0002.tmp.vcf.gz
  tabix --preset vcf ${isec_dir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${isec_dir}/0003.annot.vcf.gz \
    --columns +FORMAT,Strelka2FILTER \
    --output-type z \
    --output ${isec_dir}/0002.annot.vcf.gz \
    ${isec_dir}/0002.tmp.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0001.vcf.gz \
    --mark-sites +Strelka2 \
    --output-type z \
    --output ${isec_dir}/0001.annot.vcf.gz \
    ${isec_dir}/0001.vcf.gz

  tabix --preset vcf ${isec_dir}/0000.annot.vcf.gz
  tabix --preset vcf ${isec_dir}/0001.annot.vcf.gz
  tabix --preset vcf ${isec_dir}/0002.annot.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ${isec_dir}/0000.annot.vcf.gz \
    ${isec_dir}/0001.annot.vcf.gz \
    ${isec_dir}/0002.annot.vcf.gz | \
  bcftools sort | \
  bcftools annotate \
    --header-lines vcf.rm.header \
    --annotations ${repeatMasker} \
    --columns CHROM,FROM,TO,RepeatMasker | \
  bcftools annotate \
    --header-lines vcf.map.header \
    --annotations ${mapabilityBlacklist} \
    --columns CHROM,FROM,TO,EncodeDacMapability \
    --output-type z \
    --output ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.vcf.gz

  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    --output-type z \
    --output ${idTumor}.union.gnomad.vcf.gz \
    ${idTumor}.union.vcf.gz

  tabix ${idTumor}.union.gnomad.vcf.gz

  bcftools annotate \
    --header-lines vcf.pon.header \
    --annotations ${pon} \
    --columns PoN:=AC_Het \
    ${idTumor}.union.gnomad.vcf.gz | \
    vt annotate_indels \
    -r ${genomeFile} \
    -o ${idTumor}.union.annot.vcf -
  
  filter-vcf.py ${idTumor}.union.annot.vcf

  mv ${idTumor}.union.annot.filter.vcf ${idTumor}_vs_${idNormal}.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}_vs_${idNormal}.pass.vcf \
    ${idTumor}_vs_${idNormal}.vcf
  """
}

process SomaticAnnotateMaf {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutations", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutput 
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.vepCache,
      referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFile

  when: "mutect2" in tools && "strelka2" in tools

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}.somatic"
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 10 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax,Ref_Tri \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0

  python /usr/bin/oncokb_annotator/MafAnnotator.py \
    -i ${outputPrefix}.raw.maf \
    -o ${outputPrefix}.raw.oncokb.maf

  filter-somatic-maf.R ${outputPrefix}.raw.oncokb.maf ${outputPrefix}
  """
}


// --- Run FACETS
(bamFilesForSnpPileup, bamFiles) = bamFiles.into(2)
 
process DoSnpPileup {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/facets", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamFilesForSnpPileup
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set assay, target, idTumor, idNormal, file("${outfile}") into SnpPileup

  when: 'facets' in tools 

  script:
  outfile = idTumor + "_" + idNormal + ".snp_pileup.dat.gz"
  """
  snp-pileup \
    --count-orphans \
    --pseudo-snps=50 \
    --gzip \
    ${facetsVcf} \
    ${outfile} \
    ${bamTumor} ${bamNormal}
  """
}

// FACETS

process DoFacets {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/facets", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(snpPileupFile) from SnpPileup

  output:
    set idTumor, idNormal, target, file("${outputDir}/*purity.out"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.seg"), file("${outputDir}/*hisens.out"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.seg"), file("${outputDir}/*purity.CNCF.png"), file("${outputDir}/*hisens.CNCF.png") into FacetsOutput

  when: 'facets' in tools 

  script:
  tag = "${idTumor}_vs_${idNormal}"
  countsFile = "${snpPileupFile}"
  outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"
  """
  mkdir ${outputDir}
  /usr/bin/facets-suite/doFacets.R \
    --cval ${params.facets.cval} \
    --snp_nbhd ${params.facets.snp_nbhd} \
    --ndepth ${params.facets.ndepth} \
    --min_nhet ${params.facets.min_nhet} \
    --purity_cval ${params.facets.purity_cval} \
    --purity_snp_nbhd ${params.facets.purity_snp_nbhd} \
    --purity_ndepth ${params.facets.purity_ndepth} \
    --purity_min_nhet ${params.facets.purity_min_nhet} \
    --genome ${params.facets.genome} \
    --counts_file ${countsFile} \
    --TAG ${tag} \
    --directory ${outputDir} \
    --R_lib /usr/lib/R/library \
    --seed ${params.facets.seed} \
    --tumor_id ${idTumor}
  """
}

// --- Run MSIsensor
(bamsForMsiSensor, bamFiles) = bamFiles.into(2)

process RunMsiSensor {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/msisensor", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    file("${idTumor}_vs_${idNormal}.msisensor.tsv") into msiOutput 

  when: "msisensor" in tools

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}.msisensor.tsv"
  """
  msisensor msi \
    -d ${msiSensorList} \
    -t ${bamTumor} \
    -n ${bamNormal} \
    -o ${outputPrefix}
  """
}

// --- Run Polysolver 
(bamsForPolysolver, bamFiles) = bamFiles.into(2)

process RunPolysolver {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/hla", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForPolysolver

  output:
    file("${outputDir}/winners.hla.txt") into hlaOutput

  when: "polysolver" in tools
  
  script:
  outputDir = "."
  tmpDir = "${outputDir}-nf-scratch"
  """
  cp /home/polysolver/scripts/shell_call_hla_type .
  
  sed -i "171s/TMP_DIR=.*/TMP_DIR=${tmpDir}/" shell_call_hla_type 

  bash shell_call_hla_type \
  ${bamNormal} \
  Unknown \
  1 \
  hg19 \
  STDFQ \
  0 \
  ${outputDir} || echo "HLA Polysolver did not run successfully and its process has been redirected to generate this file." > ${outputDir}/winners.hla.txt 
  """
}

// --- Run Conpair

(bamsForConpair, bamFiles) = bamFiles.into(2)

process RunConpair {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/qc/conpair", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForConpair
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set file("${idNormal}.pileup"), file("${idTumor}.pileup") into conpairPileup
    set file("${idTumor}_${idNormal}_concordance.txt"), file("${idTumor}_${idNormal}_contamination.txt") into conpairOutput

  when: 'conpair' in tools

  script:
  gatkPath = "/usr/bin/GenomeAnalysisTK.jar"
  conpairPath = "/usr/bin/conpair"

  // These marker files are in the conpair container
  markersBed = ""
  markersTxt = ""

  if(params.genome == "GRCh37") {
    markersBed = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed"
    markersTxt = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
  }
  else {
    markersBed = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
    markersTxt = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
  }

  mem = 0
  if(params.mem_per_core) {
    mem = task.memory.toString().split(" ")[0].toInteger() - 1
  }
  else {
    mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
  }

  javaMem = "${mem}g"

  """
  # Make pileup files
  ${conpairPath}/scripts/run_gatk_pileup_for_sample.py \
    --gatk=${gatkPath} \
    --bam=${bamTumor} \
    --markers=${markersBed} \
    --reference=${genomeFile} \
    --xmx_java=${javaMem} \
    --outfile=${idTumor}.pileup

  ${conpairPath}/scripts/run_gatk_pileup_for_sample.py \
    --gatk=${gatkPath} \
    --bam=${bamNormal} \
    --markers=${markersBed} \
    --reference=${genomeFile} \
    --xmx_java=${javaMem} \
    --outfile=${idNormal}.pileup

  # Make pairing file
  echo "${idNormal}\t${idTumor}" > pairing.txt

  # Verify concordance
  ${conpairPath}/scripts/verify_concordances.py \
    --tumor_pileup=${idTumor}.pileup \
    --normal_pileup=${idNormal}.pileup \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --normal_homozygous_markers_only \
    --outpre=${idTumor}_${idNormal}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${idTumor}.pileup \
    --normal_pileup=${idNormal}.pileup \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${idTumor}_${idNormal}
  """
}



// Run LOHHLA

(bamsForLOHHLA, bamFiles) = bamFiles.into(2)

// Channel currently in order [ assay, target, tumorID, normalID, tumorBam, normalBam, tumorBai, normalBai ]

// Re-order bamsForLOHHLA into idTumor, idNormal, and target, i.e. 
// [ tumorID, normalID, target, tumorBam, normalBam, tumorBai, normalBai ]

bamsForLOHHLA = bamsForLOHHLA.map{ 
  item -> 
    def assay = item[0]
    def target = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def tumorBam = item[4]
    def normalBam = item[5]
    def tumorBai = item[6]
    def normalBai = item[7]

    return [ idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai ]
  }

// Polysolver channel currently in order []
// [ idTumor, idNormal, target, winners.hla.txt ]

// FACETS channel in order
// [ idTumor, idNormal, target, file("${outputDir}/*purity.Rdata"), file("${outputDir}/*.*") ]
// [idTumor, idNormal, target, *purity.out, *purity.cncf.txt, *purity.Rdata, purity.seg, hisens.out, hisens.cncf.txt, hisens.Rdata, hisens.seg into FacetsOutput

(facetsForLOHHLA, FacetsforMafAnno, FacetsOutput) = FacetsOutput.into(3)


facetsForLOHHLA = facetsForLOHHLA.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def purity_out = item[3]
    def purity_cncf = item[4]
    def purity_rdata = item[5]
    def purity_seg = item[6]
    def hisens_out = item[7]
    def hisens_cncf = item[8]
    def hisens_rdata = item[9]
    def hisens_seg = item[10]
    def purityCNCF_png = item[11]
    def hisensCNCF_png = item[12]

    return [ idTumor, idNormal, target, purity_out ]
  }


(hlaOutputForLOHHLA, hlaOutput) = hlaOutput.into(2)

// *purity.out from FACETS, winners.hla.txt from POLYSOLVER, with the above

//apply *.groupTuple(by: [0,1,2]) in order to group the channel by idTumor, idNormal, and target

mergedChannelLOHHLA = bamsForLOHHLA.combine(hlaOutputForLOHHLA, by: [0,1,2]).combine(facetsForLOHHLA, by: [0,1,2]).unique()


process RunLOHHLA {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/lohhla", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("winners.hla.txt"), file("*_purity.out") from mergedChannelLOHHLA
    set file(hlaFasta), file(hlaDat) from Channel.value([ referenceMap.hlaFasta, referenceMap.hlaDat ])

  output:
    file("*") into lohhlaOutput

  when: "lohhla" in tools && "polysolver" in tools && "facets" in tools && runSomatic

    // NOTE: --cleanUp in LOHHLAscript.R by default set to FALSE

  script:
    """
    cat winners.hla.txt | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
    
    PURITY=\$(grep Purity *_purity.out | grep -oP "[0-9\\.]+")
    PLOIDY=\$(grep Ploidy *_purity.out | grep -oP "[0-9\\.]+")
    cat <(echo -e "tumorPurity\ttumorPloidy") <(echo -e "\$PURITY\t\$PLOIDY") > tumor_purity_ploidy.txt

    Rscript /lohhla/LOHHLAscript.R \
        --patientId ${idTumor}_vs_{idNormal} \
        --normalBAMfile ${bamNormal} \
        --tumorBAMfile ${bamTumor} \
        --HLAfastaLoc ${hlaFasta} \
        --HLAexonLoc ${hlaDat} \
        --CopyNumLoc tumor_purity_ploidy.txt \
        --hlaPath massaged.winners.hla.txt \
        --gatkDir /picard-tools \
        --novoDir /opt/conda/bin
    """
}



// --- Run Mutational Signatures, github.com/mskcc/mutation-signatures, original Alexandrov et al 2013

(mafFileForMafAnno, mafFileForMutSig, mafFile) = mafFile.into(3)

process RunMutationSignatures {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutation_signatures", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig

  output:
    file("${idTumor}_vs_${idNormal}.mutsig.txt") into mutSigOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools && "mutsig" in tools && runSomatic

  script:
  """
  python /mutation-signatures/main.py \
    /mutation-signatures/Stratton_signatures30.txt \
    ${idTumor}_vs_${idNormal}.maf \
    ${idTumor}_vs_${idNormal}.mutsig.txt
  """
}


//Formatting the channel to be: idTumor, idNormal, target, purity_rdata

FacetsforMafAnno = FacetsforMafAnno.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def purity_out = item[3]
    def purity_cncf = item[4]
    def purity_rdata = item[5]
    def purity_seg = item[6]
    def hisens_out = item[7]
    def hisens_cncf = item[8]
    def hisens_rdata = item[9]
    def hisens_seg = item[10]
    def purityCNCF_png = item[11]
    def hisensCNCF_png = item[12]
    
    return [ idTumor, idNormal, target, purity_rdata ]
  }


//Formatting the channel to be grouped by idTumor, idNormal, and target

// FacetsOutput = FacetsOutput.groupTuple(by: [0,1,2])

mafFileForMafAnno = mafFileForMafAnno.groupTuple(by: [0,1,2])

FacetsMafFileCombine = FacetsforMafAnno.combine(mafFileForMafAnno, by: [0,1,2]).unique()


process DoMafAnno {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets_maf", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(facetsFiles), file(maf) from FacetsMafFileCombine

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.facets.maf") into MafAnnoOutput

  when: 'facets' in tools && "mutect2" in tools && "manta" in tools && "strelka2" in tools && runSomatic

  script:
  mapFile = "${idTumor}_${idNormal}.map"
  """
  echo "Tumor_Sample_Barcode\tRdata_filename" > ${mapFile}
  echo "${idTumor}\t${purity_rdata.fileName}" >> ${mapFile}

  /usr/bin/facets-suite/mafAnno.R -f ${mapFile} -m ${maf} -o ${idTumor}_vs_${idNormal}.facets.maf
  """
}

(mafFileForNeoantigen, mafFile) = mafFile.into(2)
mafFileForNeoantigen = mafFileForNeoantigen.groupTuple(by: [0,1,2])

hlaOutput = hlaOutput.combine(mafFileForNeoantigen, by: [0,1,2]).unique()

process RunNeoantigen {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/neoantigen", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(polysolverFile), file(mafFile) from hlaOutput
    set file(neoantigenCDNA), file(neoantigenCDS) from Channel.value([
      referenceMap.neoantigenCDNA,
      referenceMap.neoantigenCDS
    ]) // These reference files are in the neoantigen-docker.config

  output:
    set idTumor, idNormal, target, file("${outputDir}/*") into neoantigenOut

  when: "neoantigen" in tools

  // must set full path to tmp directories for netMHC and netMHCpan to work;
  // for some reason doesn't work with /scratch, so putting them in the process workspace
  script:
  outputDir = "neoantigen"
  tmpDir = "${outputDir}-tmp"
  tmpDirFullPath = "\$PWD/${tmpDir}/"
  """
  export TMPDIR=${tmpDirFullPath}
  mkdir -p ${tmpDir}
  chmod 777 ${tmpDir}

  python /usr/local/bin/neoantigen/neoantigen.py \
    --config_file /usr/local/bin/neoantigen/neoantigen-docker.config \
    --sample_id ${idTumor} \
    --hla_file ${polysolverFile} \
    --maf_file ${mafFile} \
    --output_dir ${outputDir}
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
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex"),
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
    // gnomAD resources
    result_array << ['gnomadWesVcf' : checkParamReturnFile("gnomadWesVcf")]
    result_array << ['gnomadWesVcfIndex' : checkParamReturnFile("gnomadWesVcfIndex")]
    // HLA FASTA and *dat for LOHHLA 
    result_array << ['hlaFasta' : checkParamReturnFile("hlaFasta")] 
    result_array << ['hlaDat' : checkParamReturnFile("hlaDat")] 
    // files for neoantigen & NetMHC
    result_array << ['neoantigenCDNA' : checkParamReturnFile("neoantigenCDNA")]
    result_array << ['neoantigenCDS' : checkParamReturnFile("neoantigenCDS")]
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
