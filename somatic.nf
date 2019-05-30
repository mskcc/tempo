#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - SomaticDellyCall
 - SomaticDellyFilter
 - CreateIntervalBeds
 - RunMutect2
 - RunMutect2Filter
 - SomaticCombineMutect2VCF
 - SomaticRunManta
 - SomaticRunStrelka
 - SomaticCombineChannel
 - SomaticRunBCFToolsFilterNorm
 - SomaticRunBCFToolsMerge
 - SomaticRunVCF2MAF
 - SomaticDoSNPPileup
 - DoFacets
 - RunMsiSensor
 - RunConpair
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
    set idTumor, idNormal, target, svType, file("${idTumor}_vs_${idNormal}_${svType}.bcf"), file("${idTumor}_vs_${idNormal}_${svType}.bcf.csi") into dellyCallOutput

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

process SomaticDellyFilter {
  tag {idTumor + "_vs_" + idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/delly", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, svType, file(dellyBcf), file(dellyBcfIndex) from dellyCallOutput

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}_${svType}.filter.bcf") into dellyFilterOutput

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

//Associating interval_list files with BAM files, putting them into one channel
agilentIList = agilentIntervals.map{ n -> [ n, "agilent" ] }
idtIList = idtIntervals.map{ n -> [ n, "idt" ] }
wgsIList = wgsIntervals.map{ n -> [ n, "wgs" ] }

( aBamList, iBamList, wBamList ) = bamsForMutect2Intervals.into(3)

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

  when: 'mutect2' in tools

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
    set file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervals) from Channel.value([
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

mutectStrelkaChannel = mutect2CombinedVcfOutput.combine( strelkaOutputMerged, by: [0,1,2] ).unique()

process SomaticCombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutect2combinedVCF), file(mutect2combinedVCFIndex), file(strelkaVCF), file(strelkaVCFIndex) from mutectStrelkaChannel
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile
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
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf,
      referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf,
      referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.pass.vcf") into vcfMergedOutput

  when: 'manta' in tools && 'strelka2' in tools && 'mutect2' in tools && runSomatic

  script:
  isec_dir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWgsVcf
  if (target != 'wgs') {
    pon = exomePoN
    gnomad =gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">\n##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">\n##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" > vcf.header
  echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
  echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header
  echo -e '##INFO=<ID=PoN,Number=1,Type=Integer,Description="Count in panel of normals">' > vcf.pon.header

  bcftools isec \
    --output-type z \
    --prefix ${isec_dir} \
    ${mutect2combinedVCF} ${strelkaVCF}

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
  
  filter-script.py ${idTumor}.union.annot.vcf

  mv ${idTumor}.union.annot.filter.vcf ${idTumor}_vs_${idNormal}.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}_vs_${idNormal}.pass.vcf \
    ${idTumor}_vs_${idNormal}.vcf
  """
}

process SomaticRunVcf2Maf {
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
    set idTumor, idNormal, target, file("*.maf") into mafFile

  // when: "mutect2" in tools && "manta" in tools && "strelka2" in tools && runSomatic 

  script:
  outfile="${vcfMerged}".replaceFirst(".pass.vcf", ".maf")
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /opt/vep/src/ensembl-vep \
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
    --output-maf ${outfile} \
    --filter-vcf 0
  """
}

// --- Run FACETS
(bamFilesForSnpPileup, bamFiles) = bamFiles.into(2)
 
process DoSnpPileup {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamFilesForSnpPileup
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

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(snpPileupFile) from SnpPileup

  output:
    set idTumor, idNormal, target, file("*purity.Rdata"), file("*.*") into FacetsOutput

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

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/hla_polysolver", mode: params.publishDirMode

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
    set file("${idTumor}.${idNormal}_concordance.txt"), file("${idTumor}.${idNormal}_contamination.txt") into conpairOutput

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
    --outpre=${idTumor}.${idNormal}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${idTumor}.pileup \
    --normal_pileup=${idNormal}.pileup \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${idTumor}.${idNormal}
  """
}

( mafFileForMafAnno, mafFileForMutSig ) = mafFile.into(2)

process RunMutationSignatures {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutation_signatures", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    file("*.maf") into mutSigOutput

  when: "mutect2" in tools && "manta" in tools && "strelka2" in tools && "mutsig" in tools

  script:
  """
  python /mutation-signatures/make_trinuc_maf.py \
    ${genomeFile} \
    ${maf} \
    ${idTumor}_${idNormal}.trinuc.maf

  python /mutation-signatures/main.py \
    /mutation-signatures/Stratton_signatures30.txt \
    ${idTumor}_${idNormal}.trinuc.maf \
    ${idTumor}_${idNormal}.mutsig.maf 
  """
}

FacetsOutput = FacetsOutput.groupTuple(by: [0,1,2])

mafFileForMafAnno = mafFileForMafAnno.groupTuple(by: [0,1,2])

FacetsMafFileCombine = FacetsOutput.combine(mafFileForMafAnno, by: [0,1,2]).unique()

process DoMafAnno {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets_maf", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(purityRdata), file(facetsFiles), file(maf) from FacetsMafFileCombine

  output:
    set idTumor, idNormal, target, file("${idTumor}_${idNormal}.maf") into MafAnnoOutput

  when: 'facets' in tools && "mutect2" in tools && "manta" in tools && "strelka2" in tools 

  script:
  mapFile = "${idTumor}_${idNormal}.map"
  """
  echo "Tumor_Sample_Barcode\tRdata_filename" > ${mapFile}
  echo "${idTumor}\t${purityRdata.fileName}" >> ${mapFile}

  /usr/bin/facets-suite/mafAnno.R -f ${mapFile} -m ${maf} -o ${idTumor}_${idNormal}.maf  
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
