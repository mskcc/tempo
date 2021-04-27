referenceMap = defineReferenceMap()
outDir     = file(params.outDir).toAbsolutePath()
mappingFile = file(params.bamMapping, checkIfExists: true)

TempoUtils.extractBAM(mappingFile, params.assayType)
	.into{bams4MutectGermline; bams4MantaGermline; bamsForStrelkaGermline}

process GermlineRunManta {
	tag {idNormal}

	input:
	set idNormal, target, file(bamNormal), file(baiNormal) from bams4MantaGermline
	set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

	output:
	set idNormal, target, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaOutput
	
	script:
	outputPrefix = "${idNormal}"
	options = ""
	if (params.assayType == "exome") options = "--exome"
	"""
	configManta.py \\
	  ${options} \\
	  --callRegions ${svCallingIncludeRegions} \\
	  --reference ${genomeFile} \\
	  --bam ${bamNormal} \\
	  --runDir Manta
	python Manta/runWorkflow.py \\
	  --mode local \\
	  --jobs ${task.cpus}

	mv Manta/results/variants/candidateSmallIndels.vcf.gz \\
	  Manta_${outputPrefix}.candidateSmallIndels.vcf.gz
	mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
	  Manta_${outputPrefix}.candidateSmallIndels.vcf.gz.tbi
	"""
}

bamsForStrelkaGermline.join(mantaOutput, by:[0,1])
	.set{bamsForStrelkaGermline}

process GermlineRunStrelka2 {
	tag {idNormal}

	publishDir "${outDir}/pon/${idNormal}/strelka2", mode: params.publishDirMode

	input:
	set idNormal, target, file(bamNormal), file(baiNormal), file(mantaCSI), file(mantaCSI_idx) from bamsForStrelkaGermline
	set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(idtTargetsIndex) from Channel.value([referenceMap.idtTargets, referenceMap.idtTargetsIndex])
    set file(idtv2Targets), file(idtv2TargetsIndex) from Channel.value([referenceMap.idtv2Targets, referenceMap.idtv2TargetsIndex])
    set file(agilentTargets), file(agilentTargetsIndex) from Channel.value([referenceMap.agilentTargets, referenceMap.agilentTargetsIndex])
    set file(wgsIntervals), file(wgsIntervalssIndex) from Channel.value([referenceMap.wgsTargets, referenceMap.wgsTargetsIndex])

	output:
	set idNormal, target, file("${idNormal}.strelka2.vcf.gz"), file("${idNormal}.strelka2.vcf.gz.tbi") into StrelkaGermline_out
	set file("listoffiles.txt"), file("listoffiles2.txt") into outFiles

	script:
	options = ""
	intervals = wgsIntervals
	if (params.assayType == "exome") {
		options = "--exome"
		if (target == 'agilent') intervals = agilentTargets
		if (target == 'idt') intervals = idtTargets
		if (target == 'idt_v2') intervals = idtTargets
	}
"""
configureStrelkaGermlineWorkflow.py \\
    ${options} \\
    --callRegions ${intervals} \\
    --referenceFasta ${genomeFile} \\
    --bam ${bamNormal} \\
    --indelCandidates ${mantaCSI} \\
    --runDir Strelka

python Strelka/runWorkflow.py \\
    --mode local \\
    --jobs ${task.cpus}

find . -type f > listoffiles.txt

bcftools filter \\
    --include 'FORMAT/AD[0:1]>1' \\
    Strelka/results/variants/variants.vcf.gz | \\
    bcftools norm \\
    --fasta-ref ${genomeFile} \\
    -check-ref s \\
    --multiallelics -both \\
    --output-type z \\
    --output ${idNormal}.strelka2.vcf.gz

tabix --preset vcf ${idNormal}.strelka2.vcf.gz

find . -type f > listoffiles2.txt
"""
}

process CreateScatteredIntervals {
 	input:
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
      ])
    set file(idtTargets), file(idtv2Targets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(idtv2TargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.idtv2Targets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.idtv2TargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
      ])
  
	output:
    set file("agilent*.interval_list"), val("agilent"), val("agilent") into agilentIList
    set file("idt*.interval_list"), val("idt"), val("idt") into idtIList
    set file("wgs*.interval_list"), val("wgs"), val("wgs") into wgsIList
    set file("idt*.interval_list"), val("idt_v2"), val("idt_v2") into idtv2IList

	script:
	scatterCount = params.scatterCount
	"""
	gatk SplitIntervals \\
	  --reference ${genomeFile} \\
	  --intervals ${agilentTargets} \\
	  --scatter-count ${scatterCount} \\
	  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
	  --output agilent
	for i in agilent/*.interval_list;
	do
		BASENAME=`basename \$i`
		mv \$i agilent-\$BASENAME
	done
	gatk SplitIntervals \\
	  --reference ${genomeFile} \\
	  --intervals ${idtTargets} \\
	  --scatter-count ${scatterCount} \\
	  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
	  --output idt
	for i in idt/*.interval_list;
	do
		BASENAME=`basename \$i`
		mv \$i idt-\$BASENAME
	done
	gatk SplitIntervals \\
	  --reference ${genomeFile} \\
	  --intervals ${idtv2Targets} \\
	  --scatter-count ${scatterCount} \\
	  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
	  --output idt_v2
	for i in idt_v2/*.interval_list;
	do
		BASENAME=`basename \$i`
		mv \$i idt_v2-\$BASENAME
	done
	gatk SplitIntervals \\
	  --reference ${genomeFile} \\
	  --intervals ${wgsTargets} \\
	  --scatter-count ${scatterCount} \\
	  --subdivision-mode INTERVAL_SUBDIVISION \\
	  --output wgs 
	for i in wgs/*.interval_list;
	do
		BASENAME=`basename \$i`
		mv \$i wgs-\$BASENAME
	done
  """
}

agilentIList.mix(idtIList, wgsIList, idtv2IList).set{mergedIList4N}

bams4MutectGermline.combine(mergedIList4N, by:1)
.map{ item ->
    def idNormal = item[1]
    def target = item[0]
    def normalBam = item[2]
    def normalBai = item[3]
    def intervalBed = item[4]
    def key = "${idNormal}@${target}" // adding one unique key

    return [ key, idNormal, target, normalBam, normalBai, intervalBed ]
}.map{ 
    key, idNormal, target, normalBam, normalBai, intervalBed -> 
    tuple ( 
         groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
         idNormal, target, normalBam, normalBai, intervalBed
    )
}
.transpose()
.set{ mergedChannelGermline }

process GermlineRunMutect2 {
	tag {idNormal}

	input:
	set key, idNormal, target, file(bamNormal), file(baiNormal), file(intervalBed) from mergedChannelGermline 
	set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

	output:
    set id, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forGermlineMutect2Combine

	script:
	mutect2Vcf = "${idNormal}_${intervalBed.baseName}.vcf.gz"
  	prefix = "${mutect2Vcf}".replaceFirst(".vcf.gz", "")
	"""
	gatk --java-options -Xmx8g Mutect2 \\
    --reference ${genomeFile} \\
    --intervals ${intervalBed} \\
    --input ${bamNormal} \\
     --tumor ${idNormal} \\
     --output ${mutect2Vcf}

    gatk --java-options -Xmx8g FilterMutectCalls \\
    --variant ${mutect2Vcf} \\
    --stats ${prefix}.Mutect2FilteringStats.tsv \\
    --output ${prefix}.gatk-filtered.vcf.gz
	
	find . -type f > listoffiles.txt
	"""
}

forGermlineMutect2Combine.groupTuple(by:[0,1,2]).set{ forGermlineMutect2Combine }

process GermlineCombineMutect2Vcf {
  tag {idNormal}

  publishDir "${outDir}/pon/${idNormal}/mutect2", mode: params.publishDirMode

  input:
    set id, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forGermlineMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idNormal, target, file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcf4Combine

  when: "mutect2" in tools 

  script:
  //def idNormal = id.toString().split("@")[0].split("__")[1]
  //def target   = id.toString().split("@")[1]
  def outfile  = "${idNormal}.mutect2.vcf.gz"
  """
  bcftools concat \\
    --allow-overlaps \\
    ${mutect2Vcf} | \\
  bcftools sort | \\
  bcftools norm \\
    --fasta-ref ${genomeFile} \\
    --check-ref s \\
    --multiallelics -both | \\
  bcftools norm --rm-dup all | \\
  bcftools view \\
    --samples ${idNormal} \\
    --output-type z | \\
  bcftools filter \\
    --include 'FORMAT/AD[0:1]>1' | \\
    sed -e 's/ID=RU,Number=1/ID=RU,Number=A/' -e 's/ID=AD,Number=R/ID=AD,Number=./' | \\
  bcftools norm \\
    --fasta-ref ${genomeFile} \\
    --check-ref s \\
    --multiallelics -both \\
    --output-type z \\
    --output ${outfile}

  tabix --preset vcf ${outfile}
  """
}

mutect2CombinedVcf4Combine.groupTuple(by:[1])
	.combine(StrelkaGermline_out.groupTuple(by:[1]), by:[1])
	.set{ponInputs}

process generatePoN {
	tag {target}

	publishDir "${outDir}/pon/${idNormal}/pon/${target}", mode: params.publishDirMode
	
	input:
	set idNormals, target, file(mutectCalls), file(mutectCallsTbi), file(strelkaCalls), file(strelkaCallsTbi) from ponInputs

	output:
	file("pon.annot.vcf.gz") into ponOutput

	script:
	"""
	bcftools merge \
		--merge none \
		--output-type z \
		--output pon.vcf.gz \
		pon/*vcf.gz

	bcftools +fill-tags pon.vcf.gz \\
		--output-type z \\
		--output pon.annot.vcf.gz \\
		--tags AC

	tabix --preset vcf pon.annot.vcf.gz
	"""

}

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  if(params."${item}" == null){println "${item} is not found in reference map"; exit 1}
  if(file(params."${item}", checkIfExists: false) == []){println "${item} is not found; glob pattern produces empty list"; exit 1}
  return file(params."${item}", checkIfExists: true)
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  result_array = [
    // genome reference dictionary
    'genomeDict' : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile' : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex' : checkParamReturnFile("genomeIndex"),
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'svCallingExcludeRegions' : checkParamReturnFile("svCallingExcludeRegions"),
    'svCallingIncludeRegions' : checkParamReturnFile("svCallingIncludeRegions"),
    'svCallingIncludeRegionsIndex' : checkParamReturnFile("svCallingIncludeRegionsIndex"),
    // Target and Bait BED files
    'idtTargets' : checkParamReturnFile("idtTargets"),
    //'idtTargetsUnzipped' : checkParamReturnFile("idtTargetsUnzipped"),
    'idtTargetsIndex' : checkParamReturnFile("idtTargetsIndex"),
    'idtTargetsList' : checkParamReturnFile("idtTargetsList"),  
    'idtBaitsList' : checkParamReturnFile("idtBaitsList"), 
    'idtv2Targets' : checkParamReturnFile("idtv2Targets"),
    'idtv2TargetsIndex' : checkParamReturnFile("idtv2TargetsIndex"),
    'idtv2TargetsList' : checkParamReturnFile("idtv2TargetsList"),  
    'idtv2BaitsList' : checkParamReturnFile("idtv2BaitsList"), 
    'agilentTargets' : checkParamReturnFile("agilentTargets"),
    //'agilentTargetsUnzipped' : checkParamReturnFile("agilentTargetsUnzipped"),
    'agilentTargetsIndex' : checkParamReturnFile("agilentTargetsIndex"),
    'agilentTargetsList' : checkParamReturnFile("agilentTargetsList"),  
    'agilentBaitsList' : checkParamReturnFile("agilentBaitsList"), 
    'wgsTargets' : checkParamReturnFile("wgsTargets"),
    //'wgsTargetsUnzipped' : checkParamReturnFile("wgsTargetsUnzipped"),
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex")
  ]

  if (workflow.profile != "test") {
    result_array << ['vepCache' : checkParamReturnFile("vepCache")]
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
    result_array << ['idtCodingBed' : checkParamReturnFile("idtCodingBed")]
    result_array << ['idtv2CodingBed' : checkParamReturnFile("idtv2CodingBed")]
    result_array << ['agilentCodingBed' : checkParamReturnFile("agilentCodingBed")]    
    result_array << ['wgsCodingBed' : checkParamReturnFile("wgsCodingBed")]  
  }
  return result_array
}
