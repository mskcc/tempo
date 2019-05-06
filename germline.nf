#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - DellyCall
 - DellyFilter
 - CreateIntervalBeds
 - RunHaplotypecaller
 - RunManta
 - RunStrelka
 - RunBcfToolsFilterNorm
 - RunBcfToolsMerge
 - RunVcf2Maf
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
  tag {idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/delly"

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, svType, file("${idNormal}_${svType}.bcf"), file("${idNormal}_${svType}.bcf.csi") into dellyCallOutput

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

process DellyFilter {
  tag {idNormal + '_' + svType}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/delly"

  input:
    set idTumor, idNormal, svType, file(dellyBcf), file(dellyBcfIndex) from dellyCallOutput


  output:
    file("*.filter.bcf") into dellyFilterOutput
    file("*.filter.bcf.csi") into dellyFilterIndexedOutput

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
(sampleIdsForIntervalBeds, bamFiles) = bamFiles.into(2)

process CreateScatteredIntervals {
  tag {idNormal}

  // publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/intervals"

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

  when: "haplotypecaller" in tools

  script:
  intervals = wgsIntervals
  if(params.assayType == "exome") {
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

(bamsForHaplotypecaller, bamFiles) = bamFiles.into(2)
bamsForHaplotypecallerIntervals = bamsForHaplotypecaller.spread(bedIntervals)

if (params.verbose) bamsForHaplotypecallerIntervals = bamsForHaplotypecallerIntervals.view {
  "BAMs for Haplotypecaller with Intervals:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  File  : [${it[4].fileName}]"
}

process RunHaplotypecaller {
  tag {idTumor + "_vs_" + idNormal + "_" + intervalBed.baseName}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/haplotypecaller"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from bamsForHaplotypecallerIntervals
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    file("*.vcf.gz") into haplotypecallerOutput
    file("*.vcf.gz.tbi") into haplotypecallerIndexedOutput

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
  """
}

(sampleIdsForHaplotypecallerCombine, bamFiles) = bamFiles.into(2)

process CombineHaplotypecallerVcf {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/haplotypecaller"

  input:
    file(haplotypecallerVcf) from haplotypecallerOutput.collect()
    file(haplotypecallerVcfIndex) from haplotypecallerIndexedOutput.collect()
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForHaplotypecallerCombine

  output:
    file("${outfile}") into haplotypecallerCombinedVcfOutput
    file("${outfile}.tbi") into haplotypecallerCombinedVcfOutputIndex

  when: 'haplotypecaller' in tools

  outfile="${idNormal}.haplotypecaller.filtered.vcf.gz"

  script:
  """
  bcftools concat \
    --allow-overlaps \
    ${haplotypecallerVcf} | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --multiallelics -both | \
  bcftools norm --rm-dup all \
    --samples ${idNormal},${idTumor} \
    --output z \
    --output-file ${outfile}

  tabix --preset vcf ${outfile}
  """
}

// --- Run Manta
(bamsForManta, bamsForStrelka, bamFiles) = bamFiles.into(3)

process RunManta {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/manta"

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
process RunStrelka2 {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/strelka2"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelka
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
    set file("*_genome.vcf.gz"), file("*_genome.vcf.gz.tbi") into strelkaOutputGenome
    set file("*_variants.vcf.gz"), file("*_variants.vcf.gz.tbi") into strelkaOutputVariants

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

  mv Strelka/results/variants/genome.*.vcf.gz Strelka_${idNormal}_genome.vcf.gz
  mv Strelka/results/variants/genome.*.vcf.gz.tbi Strelka_${idNormal}_genome.vcf.gz.tbi
  mv Strelka/results/variants/variants.vcf.gz Strelka_${idNormal}_variants.vcf.gz
  mv Strelka/results/variants/variants.vcf.gz.tbi Strelka_${idNormal}_variants.vcf.gz.tbi
  """
}

(sampleIdsForStrelkaMerge, bamFiles) = bamFiles.into(2)

process MergeStrelka2Vcfs {
  tag {idTumor + "_vs_" + idNormal}

  input: 
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForStrelkaMerge
    set file(strelkaGenome), file(strelkaGenomeIndex) from strelkaOutputGenome
    set file(strelkaVariants), file(strelkaVariantsIndex) from strelkaOutputVariants
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set file('*.vcf.gz'), file('*.vcf.gz.tbi') into strelkaOutput

  when: 'strelka2' in tools

  script:
  prefix = "${strelkaGenome}".replaceFirst(".vcf.gz", "")
  outfile = "${prefix}.filtered.vcf.gz"
  """
  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    ${strelkaGenome} ${strelkaVariants} | \
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

( sampleIdsForCombineChannel, bamFiles ) = bamFiles.into(2)

process CombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    file(haplotypecallercombinedVCF) from haplotypecallerCombinedVcfOutput
    file(haplotypecallercombinedVCFIndex) from haplotypecallerCombinedVcfOutputIndex
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForCombineChannel
    set file(strelkaVCF), file(strelkaVCFIndex) from strelkaOutput
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker,
      referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist,
      referenceMap.mapabilityBlacklistIndex
    ])

  output:
    file("${idTumor}.union.pass.vcf") into vcfMergedOutput

  when: 'strelka2' in tools && 'haplotypecaller' in tools

  script:  
  isec_dir = "${idTumor}.isec"
  """
  echo -e "##INFO=<ID=HaplotypeCaller,Number=0,Type=Flag,Description=\"Variant was called by HaplotypeCaller\">\n##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" > vcf.header
  echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
  echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header

  bcftools isec \
    --output-type z \
    --prefix ${isec_dir} \
    ${haplotypecallercombinedVCF} ${strelkaVCF}

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0002.vcf.gz \
    --mark-sites \"+HaplotypeCaller;Strelka2\" \
    --output-type z \
    --output ${isec_dir}/0002.annot.vcf.gz \
    ${isec_dir}/0002.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isec_dir}/0000.vcf.gz \
    --mark-sites +HaplotypeCaller \
    --output-type z \
    --output ${isec_dir}/0000.annot.vcf.gz \
    ${isec_dir}/0000.vcf.gz

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
    --output-type v \
    --output ${idTumor}.union.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}.union.pass.vcf \
    ${idTumor}.union.vcf
  """
}

(sampleIdsForVcf2Maf, bamFiles) = bamFiles.into(2)

process RunVcf2Maf {
  tag { idTumor + "_" + idNormal }

  publishDir "${ params.outDir }/${idTumor}_vs_${idNormal}/germline_variants/mutations"

  input:
    file(vcfMerged) from vcfMergedOutput
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForVcf2Maf
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.vepCache,
      referenceMap.isoforms
    ])

  output:
    file("*.maf") into mafFile

  when: "strelka2" in tools && "haplotypecaller" in tools

  outfile="${vcfMerged}".replaceFirst(".vcf", ".maf")

  // both tumor-id and normal-id flags are set to idNormal since we're not processing the tumor in germline.nf
  script:
  """
  perl /opt/vcf2maf.pl \
    --maf-center MSKCC-CMO \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data ${vepCache} \
    --vep-forks 10 \
    --tumor-id ${idNormal} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idNormal} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info HaplotypeCaller,Strelka2,RepeatMasker,EncodeDacMapability \
    --custom-enst ${isoforms} \
    --output-maf ${outfile} \
    --filter-vcf 0
  """
}

// --- Process Delly and Manta VCFs 

( sampleIdsForDellyMantaMerge, bamFiles ) = bamFiles.into(2)

process MergeDellyAndManta {
  tag { idNormal }

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/vcf_merged_output"

  input:
    file(dellyFilterData) from dellyFilterOutput.collect()
    file(mantaData) from mantaOutput.collect()
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForDellyMantaMerge

  output:
    file("*.merge.vcf.gz") into vcfDellyMantaMergedOutput

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
    --output ${idNormal}.delly.manta.merge.vcf.gz \
    *.vcf.gz
  """
}

( sampleIdsForBcfToolsFilter, bamFiles ) = bamFiles.into(2)

process RunBcfToolsFilterOnDellyManta {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/vcf_output"

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from sampleIdsForBcfToolsFilter
    file(vcf) from vcfDellyMantaMergedOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    file("*filtered.vcf.gz") into vcfFilterDellyMantaOutput

  when: "manta" in tools && "delly" in tools

  outfile = "${vcf}".replaceFirst('vcf.gz', 'filtered.vcf.gz')

  script:
  """
  tabix --preset vcf ${vcf}

  bcftools filter \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --output-type z \
    --output ${outfile} \
    ${vcf} 
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
  }
  return result_array
}

def extractBamFiles(tsvFile) {
  // Channeling the TSV file containing FASTQ.
  // Format is: "assay targets idTumor idNormal bamTumor bamNormal baiTumor baiNormal"
  Channel.from(tsvFile)
  .splitCsv(sep: '\t')
  .map { row ->
    checkNumberOfItem(row, 8)
    def assay = row[0]
    def target = row[1]
    def idTumor = row[2]
    def idNormal = row[3]
    def bamTumor = returnFile(row[4])
    def bamNormal = returnFile(row[5])
    def baiTumor = returnFile(row[6])
    def baiNormal = returnFile(row[7])

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
