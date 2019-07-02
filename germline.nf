#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - GermlineDellyCall
 - GermlineDellyFilter
 - CreateIntervalBeds
 - GermlineRunHaplotypecaller
 - GermlineRunManta
 - GermlineRunStrelka
 - GermlineRunBcfToolsFilterNorm
 - GermlineRunBcfToolsMerge
 - GermlineRunVcf2Maf
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
(bamsForDellyGermline, bamFiles) = bamFiles.into(2)

process GermlineDellyCall {
  tag {idNormal + '_' + svType}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/delly"

  input:
    each svType from Channel.from("DUP", "BND", "DEL", "INS", "INV")
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDellyGermline
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, svType, file("${idNormal}_${svType}.bcf"), file("${idNormal}_${svType}.bcf.csi") into dellyCallOutputGermline

  when: 'delly' in tools

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idNormal}_${svType}.bcf \
    ${bamNormal}
  """
}

process GermlineDellyFilter {
  tag {idNormal + '_' + svType}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/delly"

  input:
    set idTumor, idNormal, target, svType, file(dellyBcf), file(dellyBcfIndex) from dellyCallOutputGermline


  output:
    set idTumor, idNormal, target, file("${idNormal}_${svType}.filter.bcf") into dellyFilterOutputGermline

  when: 'delly' in tools

  outfile="${dellyBcf}".replaceFirst(".bcf",".filter.bcf")

  script:
  """
  delly filter \
    --filter germline \
    --outfile ${outfile} \
    ${dellyBcf}
  """
}

// --- Run Haplotypecaller
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

  when: "haplotypecaller" in tools

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

(bamsForHaplotypecaller, bamFiles) = bamFiles.into(2)

//Associating interval_list files with BAM files, putting them into one channel
agilentIList = agilentIntervals.map{ n -> [ n, "agilent" ] }
idtIList = idtIntervals.map{ n -> [ n, "idt" ] }
wgsIList = wgsIntervals.map{ n -> [ n, "wgs" ] }

(aBamList, iBamList, wBamList) = bamsForHaplotypecaller.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 1).unique() 
bMergedChannel = iBamList.combine(idtIList, by: 1).unique() 
wMergedChannel = wBamList.combine(wgsIList, by: 1).unique() 

mergedChannelGermline = aMergedChannel.concat( bMergedChannel, wMergedChannel)

if (params.verbose) bamsForHaplotypecallerIntervals = bamsForHaplotypecallerIntervals.view {
  "BAMs for Haplotypecaller with Intervals:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[4].fileName}]"
}

process GermlineRunHaplotypecaller {
  tag {idNormal + "_" + intervalBed.baseName}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/haplotypecaller"

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz"),
    file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz.tbi"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz.tbi") into haplotypecallerOutput mode flatten

  when: 'haplotypecaller' in tools

  script:
  """
  # Xmx hard-coded for now due to lsf bug
  # Wrong intervals set here
  gatk --java-options -Xmx8g \
    HaplotypeCaller \
    --reference ${genomeFile} \
    --intervals ${intervalBed} \
    --input ${bamNormal} \
    --output ${idNormal}_${intervalBed.baseName}.vcf.gz

  gatk SelectVariants \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.vcf.gz \
    --select-type-to-include SNP \
    --output ${idNormal}_${intervalBed.baseName}.snps.vcf.gz

  gatk SelectVariants \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.vcf.gz \
    --select-type-to-include INDEL \
    --output ${idNormal}_${intervalBed.baseName}.indels.vcf.gz

  gatk VariantFiltration \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.snps.vcf.gz \
    --filter-expression \"QD < 2.0\" --filter-name \"QD2\" \
    --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" \
    --filter-expression \"SOR > 3.0\" --filter-name \"SOR3\" \
    --filter-expression \"FS > 60.0\" --filter-name \"FS60\" \
    --filter-expression \"MQ < 40.0\" --filter-name \"MQ40\" \
    --filter-expression \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" \
    --filter-expression \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\" \
    --output ${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz

  gatk VariantFiltration \
    --reference ${genomeFile} \
    --variant ${idNormal}_${intervalBed.baseName}.indels.vcf.gz \
    --filter-expression \"QD < 2.0\" --filter-name \"QD2\" \
    --filter-expression \"QUAL < 30.0\" --filter-name \"QUAL30\" \
    --filter-expression \"FS > 200.0\" --filter-name \"FS200\" \
    --filter-expression \"ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSum-20\" \
    --output ${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz
  """
}

//Formatting the channel to be keyed by idTumor, idNormal, and target
haplotypecallerOutput = haplotypecallerOutput.groupTuple(by: [0,1,2])

process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/haplotypecaller"

  input:
    set idTumor, idNormal, target, file(haplotypecallerSnpVcf), file(haplotypecallerSnpVcfIndex), file(haplotypecallerIndelVcf), file(haplotypecallerIndelVcfIndex) from haplotypecallerOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools

  outfile="${idNormal}.haplotypecaller.vcf.gz"

  script:
  """
  bcftools concat \
    --allow-overlaps \
    ${haplotypecallerSnpVcf} ${haplotypecallerIndelVcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all \
    --output-type z \
    --output ${outfile}

  tabix --preset vcf ${outfile}
  """
}

// --- Run Manta
(bamsForMantaGermline, bamsForStrelkaGermline, bamFiles) = bamFiles.into(3)

process GermlineRunManta {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/manta"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMantaGermline
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions,
      referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("Manta_${idNormal}.diploidSV.vcf.gz") into mantaOutputGermline mode flatten

  when: 'manta' in tools

  // flag with --exome if exome
  script:
  options = ""
  if (params.assayType == "exome") options = "--exome"
  """
  configManta.py \
    ${options} \
    --callRegions ${svCallingIncludeRegions} \
    --reference ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Manta

  python Manta/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${idNormal}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${idNormal}.diploidSV.vcf.gz.tbi
  """
}

// --- Run Strelka2
process GermlineRunStrelka2 {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/strelka2"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelkaGermline
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
    set idTumor, idNormal, target, file("Strelka_${idNormal}_variants.vcf.gz"), file("Strelka_${idNormal}_variants.vcf.gz.tbi") into strelkaOutputGermline

  when: 'strelka2' in tools
  
  script:
  options = ""
  if (params.assayType == "exome") options = "--exome"

  intervals = wgsIntervals
  if(params.assayType == "exome") {
    if(target == 'agilent') intervals = agilentTargets
    if(target == 'idt') intervals = idtTargets
  }
  """
  configureStrelkaGermlineWorkflow.py \
    ${options} \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --bam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Strelka/results/variants/variants.vcf.gz Strelka_${idNormal}_variants.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi Strelka_${idNormal}_variants.vcf.gz.tbi
  """
}

// Join HaploTypeCaller and Strelka outputs 

hcv = haplotypecallerCombinedVcfOutput.groupTuple(by: [0,1,2])

haplotypecallerStrelkaChannel = hcv.combine(strelkaOutputGermline, by: [0,1,2]).unique()

(bamsForCombineChannel, bamFiles) = bamFiles.into(2)

bamsForCombineChannel = bamsForCombineChannel.map{
  item -> 
    def assay = item[0]
    def target = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def bamTumor = item[4]
    def bamNormal = item[5]
    def baiTumor = item[6]
    def baiNormal = item[7]
    
    return [idTumor, idNormal, target, assay, bamTumor, baiTumor]
  }

mergedChannelVcfCombine = bamsForCombineChannel.combine(haplotypecallerStrelkaChannel, by: [0,1,2]).unique()

process GermlineCombineChannel {
  tag {idNormal}

  input:
    set idTumor, idNormal, target, assay, file(bamTumor), file(baiTumor), file(haplotypecallercombinedVcf), file(haplotypecallercombinedVcfIndex), file(strelkaVcf), file(strelkaVcfIndex) from mergedChannelVcfCombine
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker,
      referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist,
      referenceMap.mapabilityBlacklistIndex
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf,
      referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf,
      referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.germline.vcf") into vcfMergedOutputGermline

  when: 'strelka2' in tools && 'haplotypecaller' in tools

  script:  
  isec_dir = "${idNormal}.isec"
  if (target == 'wgs') {
    gnomad = gnomadWgsVcf
  }
  else {
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=HaplotypeCaller,Number=0,Type=Flag,Description=\"Variant was called by HaplotypeCaller\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
  echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header

  bcftools isec \
    --output-type z \
    --prefix ${isec_dir} \
    ${haplotypecallercombinedVcf} ${strelkaVcf}

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
    --mark-sites +HaplotypeCaller \
    --output-type z \
    --output ${isec_dir}/0000.annot.vcf.gz \
    ${isec_dir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0002.vcf.gz \
    --mark-sites \"+HaplotypeCaller;Strelka2\" \
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
    --output ${idNormal}.union.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.union.pass.vcf.gz \
    ${idNormal}.union.vcf.gz

  tabix --preset vcf ${idNormal}.union.pass.vcf.gz

  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    ${idNormal}.union.pass.vcf.gz | \
  bcftools filter \
    --exclude \"non_cancer_AF_popmax>0.02\" \
    --output-type v \
    --output ${idNormal}.union.gnomad.vcf 

  GetBaseCountsMultiSample \
    --fasta ${genomeFile} \
    --bam ${idTumor}:${bamTumor} \
    --vcf ${idNormal}.union.gnomad.vcf \
    --output ${idTumor}.genotyped.vcf

  bgzip ${idNormal}.union.gnomad.vcf
  bgzip ${idTumor}.genotyped.vcf
  tabix --preset vcf ${idNormal}.union.gnomad.vcf.gz
  tabix --preset vcf ${idTumor}.genotyped.vcf.gz

  bcftools merge \
    --output ${idTumor}_vs_${idNormal}.germline.vcf \
    --output-type v \
    ${idNormal}.union.gnomad.vcf.gz \
    ${idTumor}.genotyped.vcf.gz
  """
}

process GermlineAnnotateMaf {
  tag {idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/mutations"

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutputGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.vepCache,
      referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFileGermline

  when: "strelka2" in tools && "haplotypecaller" in tools

  // both tumor-id and normal-id flags are set to idNormal since we're not processing the tumor in germline.nf
  script:
  outputPrefix = "${idTumor}_vs_${idNormal}.germline"
  if (target == 'wgs') {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,AC,AF,AC_nfe_seu,AF_nfe_seu,AC_afr,AF_afr,AC_nfe_onf,AF_nfe_onf,AC_amr,AF_amr,AC_eas,AF_eas,AC_nfe_nwe,AF_nfe_nwe,AC_nfe_est,AF_nfe_est,AC_nfe,AF_nfe,AC_fin,AF_fin,AC_asj,AF_asj,AC_oth,AF_oth,AC_popmax,AN_popmax,AF_popmax"
  }
  else {
    infoCols = "MuTect2,Strelka2,Strelka2FILTER,RepeatMasker,EncodeDacMapability,PoN,Ref_Tri,gnomAD_FILTER,non_cancer_AC_nfe_onf,non_cancer_AF_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AF_nfe_seu,non_cancer_AC_eas,non_cancer_AF_eas,non_cancer_AC_asj,non_cancer_AF_asj,non_cancer_AC_afr,non_cancer_AF_afr,non_cancer_AC_amr,non_cancer_AF_amr,non_cancer_AC_nfe_nwe,non_cancer_AF_nfe_nwe,non_cancer_AC_nfe,non_cancer_AF_nfe,non_cancer_AC_nfe_swe,non_cancer_AF_nfe_swe,non_cancer_AC,non_cancer_AF,non_cancer_AC_fin,non_cancer_AF_fin,non_cancer_AC_eas_oea,non_cancer_AF_eas_oea,non_cancer_AC_raw,non_cancer_AF_raw,non_cancer_AC_sas,non_cancer_AF_sas,non_cancer_AC_eas_kor,non_cancer_AF_eas_kor,non_cancer_AC_popmax,non_cancer_AF_popmax"
  }
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /usr/bin/vep \
    --vep-data ${vepCache} \
    --vep-forks 4 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0

  filter-germline-maf.R ${outputPrefix}.raw.maf ${outputPrefix}
  """
}

// --- Process Delly and Manta VCFs 

mantaOutputGermline = mantaOutputGermline.groupTuple(by: [0,1,2])
dellyFilterOutputGermline = dellyFilterOutputGermline.groupTuple(by: [0,1,2])

dellyMantaChannelGermline = dellyFilterOutputGermline.combine(mantaOutputGermline, by: [0,1,2]).unique()

process GermlineMergeDellyAndManta {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/vcf_merged_output"

  input:
    set idTumor, idNormal, target, file(dellyBcf), file(mantaVcf) from dellyMantaChannelGermline

  output:
    set idTumor, idNormal, target, file("${idNormal}.delly.manta.unfiltered.vcf.gz") into vcfDellyMantaMergedOutputGermline

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
    --output-type z \
    --output ${idNormal}.delly.manta.unfiltered.vcf.gz \
    *.vcf.gz
  """
}

process GermlineRunBcfToolsFilterOnDellyManta {
  tag {idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/structural_variants"

  input:
    set idTumor, idNormal, target, file(vcf) from vcfDellyMantaMergedOutputGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into vcfFilterDellyMantaOutputGermline

  when: "manta" in tools && "delly" in tools

  outfile = "${vcf}".replaceFirst('.unfiltered.vcf.gz', '.vcf.gz')

  script:
  """
  tabix --preset vcf ${vcf}

  bcftools filter \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${outfile} \
    ${vcf}

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
    // gnomAD resources
    result_array << ['gnomadWesVcf' : checkParamReturnFile("gnomadWesVcf")]
    result_array << ['gnomadWesVcfIndex' : checkParamReturnFile("gnomadWesVcfIndex")]
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
