#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
Processes overview:

Alignment and QC
----------------
 - AlignReads
    --- Map paired-end FASTQs with bwa mem
    --- Sort BAM with samtools sort
    --- FASTQ QC with FastP 
 - MergeBam --- Merge BAM for the same samples from different lanes, samtools merge
 - MarkDuplicates --- Mark Duplicates with GATK4 MarkDuplicates
 - CreateRecalibrationTable --- Create Recalibration Table with GATK4 BaseRecalibrator
 - RecalibrateBam --- Recalibrate Bam with GATK4 ApplyBQSR

Somatic Analysis
----------------
 - CreateScatteredIntervals --- GATK4 SplitIntervals
 - CollectHsMetrics --- *For WES only* Calculate hybrid-selection metrics, GATK4 CollectHsMetrics
 - RunMutect2 --- somatic SNV calling, MuTect2
 - SomaticRunStrelka2 --- somatic SNV calling, Strelka2, using Manta for small InDel calling by default
 - SomaticCombineMutect2Vcf --- combine Mutect2 calls,bcftools
 - SomaticRunManta --- somatic SV calling, Manta
 - SomaticDellyCall --- somatic SV calling, Delly
 - SomaticMergeDellyAndManta --- combine Manta and Delly VCFs
 - SomaticCombineChannel --- combine and filter VCFs, bcftools
 - SomaticAnnotateMaf --- annotate MAF, vcf2maf
 - RunMsiSensor --- MSIsensor
 - DoFacets --- facets-suite: mafAnno.R, geneLevel.R, armLevel.R
 - RunPolysolver --- Polysolver
 - RunLOHHLA --- LOH in HLA
 - RunConpair --- Tumor-Normal quality/contamination
 - RunMutationSignatures --- mutational signatures
 - FacetsAnnotation --- annotate FACETS
 - RunNeoantigen --- NetMHCpan 4.0
 - MetaDataParser --- python script to parse metadata into single *tsv
 - SomaticAggregate --- collect outputs

Germline Analysis
-----------------
 - CreateScatteredIntervals --- (run once) GATK4 SplitIntervals
 - GermlineDellyCall --- germline SV calling and filtering, Delly
 - GermlineRunManta --- germline SV calling, Manta
 - GermlineMergeDellyAndManta --- merge SV calls from Delly and Manta
 - GermlineRunHaplotypecaller --- germline SNV calling, GATK4
 - GermlineCombineHaplotypecallerVcf --- concatenate VCFs of GATK4 HaplotypeCaller
 - GermlineRunStrelka2 --- germline SNV calling, Strelka2 (with InDels from Manta)
 - GermlineCombineChannel --- combined and filter germline calls, bcftools
 - GermlineAnnotateMaf--- annotate MAF, vcf2maf
 - GermlineAggregate --- collect outputs

*/

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println "ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)"
  exit 1
}

// Both mapping and pairing necessary for alignment of FASTQs
// Only bam_pairing required when using already aligned BAM files
if (params.mapping && !params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --pairing and re-run the pipeline."
  exit 1
}

if (!params.mapping && params.pairing) {
  println "ERROR: Flags --mapping and --pairing must both be provided. Please provide --mapping and re-run the pipeline."
  exit 1
}

if ((params.mapping && params.bam_pairing) || (params.pairing && params.bam_pairing)) {
  println "ERROR: Cannot use both FASTQs and BAMs as inputs. Flags --bam_pairing and --mapping/-pairing cannot be invoked together. Please provide either FASTQs or BAMs, and re-run the pipeline."
  exit 1
} 

// Validate mapping file
// Check for duplicate inputs, mixed assay types and unique lane names
if (params.mapping) {
  mappingPath = params.mapping
  
  if (mappingPath && !check_for_duplicated_rows(mappingPath)) {
    println "ERROR: Duplicated row found in mapping file. Please fix the error and re-run the pipeline."
    exit 1
  }

  if (mappingPath && !check_for_mixed_assay(mappingPath)) {
    println "ERROR: Multiple assays found in mapping file. Users can either run exomes or genomes, but not both. Please fix the error and re-run the pipeline."
    exit 1
  }

  if (mappingPath && !checkForUniqueSampleLanes(mappingPath)) {
    println "ERROR: The combination of sample ID and lane names values must be unique. Duplicate lane names for one sample cause errors. Please fix the error and re-run the pipeline."
    exit 1
  }
}

// Validate pairing file
// Check for duplicate inputs
if (params.pairing) {
  pairingPath = params.pairing

  if (!check_for_duplicated_rows(pairingPath)) {
    println "ERROR: Duplicated row found in pairing file. Please fix the error and re-run the pipeline."
    exit 1
  }
}

// Validate BAM file pairing file
// Check for duplicate inputs
if (params.bam_pairing) {
  bamPairingPath = params.bam_pairing

  if (bamPairingPath && !check_for_duplicated_rows(bamPairingPath)) {
    println "ERROR: Duplicated row found in BAM mapping file. Please fix the error and re-run the pipeline."
    exit 1
  }
}

// User-set runtime parameters
publishAll = params.publishAll
outname = params.outname
runGermline = params.germline
runSomatic = params.somatic

referenceMap = defineReferenceMap()

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

// Skip these processes if starting from aligned BAM files
if (!params.bam_pairing) {

  // Parse input FASTQ mapping and sample pairing
  fastqFiles = Channel.empty() 
  mappingFile = file(mappingPath)
  pairingFile = file(pairingPath)
  pairingTN = extractPairing(pairingFile)
  fastqFiles = extractFastq(mappingFile)

  fastqFiles.groupTuple(by:[0]).map{ key, lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets -> tuple( groupKey(key, lanes.size()), lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets)}.set{ groupedFastqs }

  groupedFastqs.into { fastPFiles; fastqFiles }
  fastPFiles = fastPFiles.transpose()
  fastqFiles = fastqFiles.transpose()

  // AlignReads - Map reads with BWA mem output SAM
  process AlignReads {
    tag {idSample + "@" + lane}   // The tag directive allows you to associate each process executions with a custom label

    publishDir "${params.outDir}/qc/fastp/${idSample}", pattern: "*.html", mode: params.publishDirMode

    input:
      set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assay, targetFile from fastqFiles
      set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

    output:
      file("*.html") into fastPResults
      set idSample, lane, file("${lane}.sorted.bam"), assay, targetFile into sortedBam

    script:
    readGroup = "@RG\\tID:${lane}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
    // Different resource requirements for AWS and LSF
    if (params.mem_per_core) { 
      mem = task.memory.toString().split(" ")[0].toInteger() - 1 
    }
    else {
      mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
    } 
    """
    set -e
    set -o pipefail
    fastp --html ${lane}.fastp.html --in1 ${fastqFile1} --in2 ${fastqFile2} 
    bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} | samtools view -Sb - > ${lane}.bam

    samtools sort -m ${mem}G -@ ${task.cpus} -o ${lane}.sorted.bam ${lane}.bam
    """
  }

  sortedBam.groupTuple().set{ groupedBam }

  groupedBam = groupedBam.map{ item -> 
    def idSample = item[0]
    def lane = item[1] //is a list
    def bam = item[2]
  
    def assayList = item[3].unique()
    def targetList = item[4].unique()

    if ((assayList.size() > 1) || (targetList.size() > 1)) {  
      println "ERROR: Multiple assays and/or targets found for ${idSample}; check inputs"
      exit 1
    }
    
    def assay = assayList[0]
    def target = targetList[0]

    [idSample, lane, bam, assay, target]
  }


  // MergeBams
  process MergeBams {
    tag {idSample}

    input:
      set idSample, lane, file(bam), assay, targetFile from groupedBam

    output:
      set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile into mergedBam

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
    """
  }

  // GATK MarkDuplicates
  process MarkDuplicates {
    tag {idSample}

    input:
      set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile from mergedBam

    output:
      set file("${idSample}.md.bam"), file("${idSample}.md.bai"), idSample, lane, assay, targetFile into duplicateMarkedBams
      set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), assay, targetFile into markDuplicatesTSV
      file ("${idSample}.bam.metrics") into markDuplicatesReport

    script:
    """
    gatk MarkDuplicates --java-options ${params.markdup_java_options} \
      --MAX_RECORDS_IN_RAM 50000 \
      --INPUT ${idSample}.merged.bam \
      --METRICS_FILE ${idSample}.bam.metrics \
      --TMP_DIR . \
      --ASSUME_SORT_ORDER coordinate \
      --CREATE_INDEX true \
      --OUTPUT ${idSample}.md.bam
    """
  }

  duplicateMarkedBams = duplicateMarkedBams.map {
      bam, bai, idSample, lane, assay, targetFile -> tag = bam.baseName.tokenize('.')[0]
      [idSample, bam, bai, assay, targetFile]
  }

  (mdBam, mdBamToJoin) = duplicateMarkedBams.into(2)

 // GATK BaseRecalibrator , CreateRecalibrationTable
 
  process CreateRecalibrationTable {
    tag {idSample}

    input:
      set idSample, file(bam), file(bai), assay, targetFile from mdBam 
      set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex) from Channel.value([
        referenceMap.genomeFile,
        referenceMap.genomeIndex,
        referenceMap.genomeDict,
        referenceMap.dbsnp,
        referenceMap.dbsnpIndex,
        referenceMap.knownIndels,
        referenceMap.knownIndelsIndex 
      ])

    output:
      set idSample, file("${idSample}.recal.table") into recalibrationTable
      set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), val("${idSample}.recal.table"), assay, targetFile into recalibrationTableTSV

    script:
    knownSites = knownIndels.collect{ "--known-sites ${it}" }.join(' ')
    """
    gatk BaseRecalibrator \
      --tmp-dir /tmp \
      --reference ${genomeFile} \
      --known-sites ${dbsnp} \
      ${knownSites} \
      --verbosity INFO \
      --input ${bam} \
      --output ${idSample}.recal.table
    """ 
  }

  recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0])

  // GATK ApplyBQSR, RecalibrateBAM
  process RecalibrateBam {
    tag {idSample}

    publishDir "${params.outDir}/bams", mode: params.publishDirMode

    input:
      set idSample, file(bam), file(bai), assay, targetFile, file(recalibrationReport) from recalibrationTable
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict 
      ])

    output:
      set idSample, file("${idSample}.bam"), file("${idSample}.bam.bai"), assay, targetFile into recalibratedBam, recalibratedBamForCollectHsMetrics, recalibratedBamForStats, recalibratedBamForOutput, recalibratedBamForOutput2
      set idSample, val("${idSample}.bam"), val("${idSample}.bai"), assay, targetFile into recalibratedBamTSV
      val(idSample) into currentSample
      file("${idSample}.bam") into currentBam
      file("${idSample}.bai") into currentBai
      val(assay) into assays
      val(targetFile) into targets

    script:
    """
    gatk ApplyBQSR \
      --reference ${genomeFile} \
      --create-output-bam-index true \
      --bqsr-recal-file ${recalibrationReport} \
      --input ${bam} \
      --output ${idSample}.bam
      
    cp -p ${idSample}.bai ${idSample}.bam.bai
    """
  }

  // set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta

  recalibratedBamForOutput.combine(pairingTN)
                          .filter { item -> // only keep combinations where sample is same as tumor pair sample
                            def idSample = item[0]
                            def sampleBam = item[1]
                            def sampleBai = item[2]
                            def assay = item[3]
                            def target = item[4]
                            def idTumor = item[5]
                            def idNormal = item[6]
                            idSample == idTumor
                          }.map { item -> // re-order the elements
                            def idSample = item[0]
                            def sampleBam = item[1]
                            def sampleBai = item[2]
                            def assay = item[3]
                            def target = item[4]
                            def idTumor = item[5]
                            def idNormal = item[6]
                            def bamTumor = sampleBam
                            def baiTumor = sampleBai

                            return [ assay, target, idTumor, idNormal, bamTumor, baiTumor ]
                          }.combine(recalibratedBamForOutput2)
                          .filter { item ->
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def baiTumor = item[5]
                            def idSample = item[6]
                            def bamNormal = item[7]
                            def baiNormal = item[8]
                            idSample == idNormal
                          }.map { item -> // re-order the elements
                            def assay = item[0]
                            def target = item[1]
                            def idTumor = item[2]
                            def idNormal = item[3]
                            def bamTumor = item[4]
                            def baiTumor = item[5]
                            def idSample = item[6]
                            def bamNormal = item[7]
                            def baiNormal = item[8]

                            return [ assay, target, idTumor, idNormal, bamTumor, bamNormal, baiTumor, baiNormal ]
                          }
                          .set { result }

  result.into { resultTsv; bamFiles }

  File file = new File(outname)
  file.newWriter().withWriter { w ->
      w << "ASSAY\tTARGET\tTUMOR_ID\tNORMAL_ID\tTUMOR_BAM\tNORMAL_BAM\tTUMOR_BAI\tNORMAL_BAI\n"
  }

  if (workflow.profile == 'awsbatch') {
      resultTsv.subscribe { Object obj ->
        file.withWriterAppend { out ->
          out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\ts3:/${obj[4]}\ts3:/${obj[5]}\ts3:/${obj[6]}\ts3:/${obj[7]}"
        }
      }
    }
  else {
    resultTsv.subscribe { Object obj ->
      file.withWriterAppend { out ->
        out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}\t${obj[4]}\t${obj[5]}\t${obj[6]}\t${obj[7]}"
      }
    }
  }
  
  // GATK CollectHsMetrics, WES only
  process CollectHsMetrics {
    tag {idSample}

    publishDir "${params.outDir}/qc/collecthsmetrics/${idSample}", mode: params.publishDirMode

    input:
      set idSample, file(bam), file(bai), assay, target from recalibratedBamForCollectHsMetrics
      set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
        referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
      ])
      set file(idtTargetsList), file(agilentTargetsList), file(idtBaitsList), file(agilentBaitsList) from Channel.value([
        referenceMap.idtTargetsList, referenceMap.agilentTargetsList, 
        referenceMap.idtBaitsList, referenceMap.agilentBaitsList      
      ])

    output:
      file("${idSample}.hs_metrics.txt") into CollectHsMetricsStats

    when: 'wes' in assay && !params.test

    script:
    baitIntervals = ""
    targetIntervals = ""
    if (target == 'agilent'){
      baitIntervals = "${agilentBaitsList}"
      targetIntervals = "${agilentTargetsList}"
    }
    if (target == 'idt'){
      bait_intervals = "${idtBaitsList}"
      target_intervals = "${idtTargetsList}"
    }
    """
    gatk CollectHsMetrics \
      --INPUT ${bam} \
      --OUTPUT ${idSample}.hs_metrics.txt \
      --REFERENCE_SEQUENCE ${genomeFile} \
      --BAIT_INTERVALS ${baitIntervals} \
      --TARGET_INTERVALS ${targetIntervals} 
    """
  }

  // Alfred, BAM QC
  ignore_read_groups = Channel.from( true , false )
  process Alfred {
    tag {idSample + "@" + "ignore_rg_" + ignore_rg }

    publishDir "${params.outDir}/qc/alfred/${idSample}", mode: params.publishDirMode
  
    input:
      each ignore_rg from ignore_read_groups
      set idSample, file(bam), file(bai), assay, target from recalibratedBam
      file(genomeFile) from Channel.value([referenceMap.genomeFile])
      set file(idtTargets), file(agilentTargets), file(idtTargetsIndex), file(agilentTargetsIndex) from Channel.value([
        referenceMap.idtTargets, referenceMap.agilentTargets,
        referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex
      ])

    output:
      set ignore_rg, idSample, file("*.tsv.gz"), file("*.tsv.gz.pdf") into bamsQCStats

    script:
    options = ""
    if (assay == "exome") {
      if (target == "agilent") options = "--bed ${agilentTargets}"
      if (target == "idt") options = "--bed ${idtTargets}"
    }
    def ignore = ignore_rg ? "--ignore" : ""
    def outfile = ignore_rg ? "${idSample}.alfred.tsv.gz" : "${idSample}.alfred.per_readgroup.tsv.gz"
    """
    alfred qc ${options} \
      --reference ${genomeFile} \
      ${ignore} \
      --outfile ${outfile} \
      ${bam} && \
      Rscript /opt/alfred/scripts/stats.R ${outfile}
    """
  }
}

/*
================================================================================
=                                SOMATIC PIPELINE                              =
================================================================================
*/

// parse --tools parameter for downstream 'when' conditionals, e.g. when: `` 'delly ' in tools
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

// Allow shorter names
if ("mutect" in tools) {
  tools.add("mutect2")
}
if ("strelka" in tools) {
  tools.add("strelka2")
}

// If using Strelka2, run Manta as well to generate candidate indels
if ("strelka2" in tools) {
  tools.add("manta")
}

// If starting with BAM files, parse BAM pairing input
if (params.bam_pairing) {
  bamFiles = Channel.empty()
  bamPairingfile = file(bamPairingPath)
  bamFiles = extractBAM(bamPairingfile)
}

// GATK SplitIntervals, CreateScatteredIntervals
process CreateScatteredIntervals {

  input:
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
      ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
      ])
  
  output:
    set file("agilent*.interval_list"), val("agilent") into agilentIList
    set file("idt*.interval_list"), val("idt") into idtIList
    set file("wgs*.interval_list"), val("wgs") into wgsIList


  when: runSomatic || runGermline

  script:
  scatterCount = params.scatterCount
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
    --intervals ${wgsTargets} \
    --scatter-count ${scatterCount} \
    --subdivision-mode INTERVAL_SUBDIVISION \
    --output wgs 

  for i in wgs/*.interval_list;
  do
    BASENAME=`basename \$i`
    mv \$i wgs-\$BASENAME
  done
  """
}


(bamsForIntervals, bamFiles) = bamFiles.into(2)


//Associating interval_list files with BAM files, putting them into one channel

(aBamList, iBamList, wBamList) = bamsForIntervals.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 1)
iMergedChannel = iBamList.combine(idtIList, by: 1)
wMergedChannel = wBamList.combine(wgsIList, by: 1)

// These will go into mutect2 and haplotypecaller

// //(mergedChannelSomatic, mergedChannelGermline) = aMergedChannel.concat( bMergedChannel, wMergedChannel).into(2) // { mergedChannelSomatic, mergedChannelGermline }

// From Nextflow Doc: However there are use cases in which each tuple has a different size depending grouping key. In this cases use the built-in function groupKey that allows you to create a special grouping key object to which it's possible to associate the group size for a given key.
// Reference: https://github.com/nextflow-io/nextflow/issues/796
// using groupKey() function to create a unique key for each TN pair, and the key also includes number of intervalBeds for each samples
// this change will allow the merging processes of each sample only wait for relative children processes from the previous step, instead of waiting for all the processes to be done
// change from .concat to .mix because .concat will wait for all the items proceeding from the first channel were emitted
(mergedChannelSomatic, mergedChannelGermline) = aMergedChannel.mix( iMergedChannel, wMergedChannel).map{
  item ->
    def key = item[2]+"_vs_"+item[3]+"@"+item[0] // adding one unique key
    def target = item[0]
    def assay = item[1]
    def idTumor = item[2]
    def idNormal = item[3]
    def tumorBam = item[4]
    def normalBam = item[5]
    def tumorBai = item[6]
    def normalBai = item[7]
    def intervalBed = item[8]

    return [ key, target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed ]
}.map{ 
    key, target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed -> 
    tuple ( 
         groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
         target, assay, idTumor, idNormal, tumorBam, normalBam, tumorBai, normalBai, intervalBed
    )
}.transpose().into(2)



// if using strelka2, one should have manta for small InDels

if('strelka2' in tools) {
  tools.add('manta')
}


// --- Run Delly
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDelly, bamFiles) = bamFiles.into(2)

process SomaticDellyCall {
  tag {idTumor + "_vs_" + idNormal + '@' + svType}

  if (publishAll) { publishDir "${params.outDir}/somatic/delly", mode: params.publishDirMode }

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDelly
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}_${svType}.filter.bcf") into dellyFilterOutput

  when: "delly" in tools && runSomatic

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

// --- Run Mutect2
process RunMutect2 {
  tag {idTumor + "_vs_" + idNormal + "@" + intervalBed.baseName}

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set id, target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelSomatic 
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set id, idTumor, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forMutect2Combine

  when: "mutect2" in tools && runSomatic

  script:
  mutect2Vcf = "${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz"
  prefix = "${mutect2Vcf}".replaceFirst(".vcf.gz", "")
  """
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
// group by groupKey(key, intervalBed.size())
forMutect2Combine = forMutect2Combine.groupTuple()

// Combine Mutect2 VCFs, bcftools
process SomaticCombineMutect2Vcf {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set id, idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex), file(mutect2Stats) from forMutect2Combine
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into mutect2CombinedVcfOutput

  when: "mutect2" in tools && runSomatic

  script:
  idTumor = id.toString().split("_vs_")[0]
  idNormal = id.toString().split("@")[0].split("_vs_")[1]
  target = id.toString().split("@")[1]
  outfile="${idTumor}_vs_${idNormal}.mutect2.filtered.vcf.gz"
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

(bamsForManta, bamsForStrelka, bamFiles) = bamFiles.into(3)

// --- Run Manta
process SomaticRunManta {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("*.vcf.gz") into mantaOutput
    set idTumor, idNormal, target, file("*.vcf.gz.tbi") into mantatbi
    set idTumor, idNormal, target, assay, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: "manta" in tools && runSomatic

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}"
  options = ""
  if (assay == "wes") options = "--exome"
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
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${outputPrefix}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${outputPrefix}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${outputPrefix}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${outputPrefix}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    Manta_${outputPrefix}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    Manta_${outputPrefix}.somaticSV.vcf.gz.tbi
  """
}

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)

dellyFilterOutput = dellyFilterOutput.groupTuple(by: [0,1,2], size: 5)

dellyMantaCombineChannel = dellyFilterOutput.combine(mantaOutput, by: [0,1,2])

// --- Process Delly and Manta VCFs 
(sampleIdsForDellyMantaMerge, bamFiles) = bamFiles.into(2)

// Merge VCFs, Delly and Manta
process SomaticMergeDellyAndManta {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(dellyBcfs), file(mantaFile) from dellyMantaCombineChannel

  output:
    file("${outputPrefix}.delly.manta.vcf.gz") into vcfDellyMantaMergedOutput

  when: tools.containsAll(["manta", "delly"]) && runSomatic

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}"
  """ 
  for f in *.bcf
  do 
    bcftools view --output-type z \$f > \${f%.bcf}.delly.vcf.gz
  done

  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done
  
  bcftools view \
    --samples ${idTumor},${idNormal} \
    --output-type z \
    --output-file Manta_${outputPrefix}.somaticSV.tmp.vcf.gz \
    Manta_${outputPrefix}.somaticSV.vcf.gz 
    
  tabix --preset vcf Manta_${outputPrefix}.somaticSV.tmp.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${outputPrefix}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz Manta_${outputPrefix}.somaticSV.tmp.vcf.gz
  
  tabix --preset vcf ${outputPrefix}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    ${outputPrefix}.delly.manta.unfiltered.vcf.gz | \
  bcftools sort \
    --output-type z \
    --output-file ${outputPrefix}.delly.manta.vcf.gz 

  tabix --preset vcf ${outputPrefix}.delly.manta.vcf.gz 
  """
}


// --- Run Strelka2

process SomaticRunStrelka2 {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/strelka2", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, assay, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsTargets),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsTargetsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])

  output:
    set idTumor, idNormal, target, file(bamNormal), file(baiNormal), file('*merged.filtered.vcf.gz'), file('*merged.filtered.vcf.gz.tbi') into strelkaOutputMerged
    set idTumor, idNormal, target, file("*indels.vcf.gz"), file("*indels.vcf.gz.tbi"), file("*snvs.vcf.gz"), file("*snvs.vcf.gz.tbi") into strelkaOutput

  when: tools.containsAll(["manta", "strelka2"]) && runSomatic

  script:
  options = ""
  intervals = wgsTargets
  if (assay == "wes") {
    options = "--exome"
    if (target == 'agilent') intervals = agilentTargets
    if (target == 'idt') intervals = idtTargets
  }
  prefix = "${idTumor}_${idNormal}_${target}.strelka.merged"
  outfile = "${prefix}.filtered.vcf.gz"
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

  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    Strelka_${idTumor}_vs_${idNormal}_somatic_indels.vcf.gz Strelka_${idTumor}_vs_${idNormal}_somatic_snvs.vcf.gz | \
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



mutectStrelkaChannel = mutect2CombinedVcfOutput.combine(strelkaOutputMerged, by: [0,1,2])

// Combined Somatic VCFs

process SomaticCombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutectCombinedVcf), file(mutectCombinedVcfIndex), file(bamNormal), file(baiNormal), file(strelkaVcf), file(strelkaVcfIndex) from mutectStrelkaChannel
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker, referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex
    ])
    set file(exomePoN), file(wgsPoN), file(exomePoNIndex), file(wgsPoNIndex) from Channel.value([
      referenceMap.exomePoN, referenceMap.wgsPoN,
      referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.pass.vcf") into vcfMergedOutput

  when: tools.containsAll(["manta", "strelka2", "mutect2"]) && runSomatic
  
  script:
  outputPrefix = "${idTumor}_vs_${idNormal}"
  isec_dir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWgsVcf
  if (target == "wgs") {
    pon = wgsPoN
    gnomad = gnomadWgsVcf
  }
  else {
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
  echo -e "##FORMAT=<ID=alt_count_raw,Number=1,Type=Integer,Description=\"Raw normal alternate allele depth\">" > vcf.ad_n.header

  # Get set differences of variant calls
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

  # Concatenate the different sets, annotate with blacklists
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

  # Add gnomAD annotation
  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    --output-type z \
    --output ${idTumor}.union.gnomad.vcf.gz \
    ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.gnomad.vcf.gz

  # Add PoN annotation and flanking sequence
  bcftools annotate \
    --header-lines vcf.pon.header \
    --annotations ${pon} \
    --columns PoN:=AC_Het \
    ${idTumor}.union.gnomad.vcf.gz | \
    vt annotate_indels \
    -r ${genomeFile} \
    -o ${idTumor}.union.annot.vcf -

  # Do custom filter annotation, then filter variants
  filter-vcf.py ${idTumor}.union.annot.vcf

  mv ${idTumor}.union.annot.filter.vcf ${outputPrefix}.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}_vs_${idNormal}.filtered.vcf \
    ${idTumor}_vs_${idNormal}.vcf

  # Add normal read count, using all reads
  GetBaseCountsMultiSample \
    --thread 4 \
    --maq 0 \
    --fasta ${genomeFile} \
    --bam ${idNormal}:${bamNormal} \
    --vcf ${idTumor}_vs_${idNormal}.filtered.vcf \
    --output ${idNormal}.genotyped.vcf 
  
  bgzip ${idTumor}_vs_${idNormal}.filtered.vcf
  bgzip ${idNormal}.genotyped.vcf
  tabix --preset vcf ${idTumor}_vs_${idNormal}.filtered.vcf.gz
  tabix --preset vcf ${idNormal}.genotyped.vcf.gz

  bcftools annotate \
    --annotations ${idNormal}.genotyped.vcf.gz \
    --header-lines vcf.ad_n.header \
    --columns FORMAT/alt_count_raw:=FORMAT/AD \
    --output-type v \
    --output ${outputPrefix}.pass.vcf \
    ${outputPrefix}.filtered.vcf.gz
  """
}

// Run vcf2maf and apply custom filters
process SomaticAnnotateMaf {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/mutations", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutput 
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.maf") into mafFile

  when: tools.containsAll(["manta", "strelka2", "mutect2"]) && runSomatic

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}.somatic"
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
    --vep-forks 10 \
    --tumor-id ${idTumor} \
    --normal-id ${idNormal} \
    --vcf-tumor-id ${idTumor} \
    --vcf-normal-id ${idNormal} \
    --input-vcf ${vcfMerged} \
    --ref-fasta ${genomeFile} \
    --retain-info ${infoCols} \
    --retain-fmt alt_count_raw \
    --custom-enst ${isoforms} \
    --output-maf ${outputPrefix}.raw.maf \
    --filter-vcf 0
    
  python /usr/bin/oncokb_annotator/MafAnnotator.py \
    -i ${outputPrefix}.raw.maf \
    -o ${outputPrefix}.raw.oncokb.maf
    
  filter-somatic-maf.R ${outputPrefix}.raw.oncokb.maf ${outputPrefix}
  """
}

(bamsForMsiSensor, bamFiles) = bamFiles.into(2)

// --- Run MSIsensor
process RunMsiSensor {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/msisensor", mode: params.publishDirMode }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.msisensor.tsv") into msiOutput 

  when: "msisensor" in tools && runSomatic

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

(msiOutputForMetaData, msiOutput) = msiOutput.into(2)

(bamFilesForSnpPileup, bamFiles) = bamFiles.into(2)

// --- Run FACETS 
process DoFacets {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/facets", mode: params.publishDirMode }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamFilesForSnpPileup
    file(facetsVcf) from Channel.value([referenceMap.facetsVcf])

  output:
    set assay, target, idTumor, idNormal, file("${outfile}") into SnpPileup
    set idTumor, idNormal, target, file("${outputDir}/*purity.out"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.seg"), file("${outputDir}/*hisens.out"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*purity.CNCF.png") into FacetsOutput
    set file("${outputDir}/*purity.seg"), file("${outputDir}/*purity.cncf.txt"), file("${outputDir}/*purity.CNCF.png"), file("${outputDir}/*purity.Rdata"), file("${outputDir}/*purity.out") into FacetsPurity
    set file("${outputDir}/*hisens.seg"), file("${outputDir}/*hisens.cncf.txt"), file("${outputDir}/*hisens.CNCF.png"), file("${outputDir}/*hisens.Rdata"), file("${outputDir}/*hisens.out") into FacetsHisens
    file("${tag}_OUT.txt") into FacetsPurityHisensOutput
    file("${outputFacetsSubdirectory}") into FacetsOutputSubdirectories

  when: "facets" in tools && runSomatic

  script:
  outfile = idTumor + "_" + idNormal + ".snp_pileup.dat.gz"
  tag = outputFacetsSubdirectory = "${idTumor}_vs_${idNormal}"
  outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"
  """
  snp-pileup \
    --count-orphans \
    --pseudo-snps=50 \
    --gzip \
    ${facetsVcf} \
    ${outfile} \
    ${bamTumor} ${bamNormal}

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
    --counts_file ${outfile} \
    --TAG ${tag} \
    --directory ${outputDir} \
    --R_lib /usr/lib/R/library \
    --seed ${params.facets.seed} \
    --tumor_id ${idTumor}

  python3 /usr/bin/facets-suite/summarize_project.py \
    -p ${tag} \
    -c ${outputDir}/*cncf.txt \
    -o ${outputDir}/*out \
    -s ${outputDir}/*seg  
  
  mkdir ${outputFacetsSubdirectory}
  cp -rf ${outfile} ${outputFacetsSubdirectory}
  cp -rf ${outputDir} ${outputFacetsSubdirectory}
  """
}

(bamsForPolysolver, bamFiles) = bamFiles.into(2)

// Run Polysolver
process RunPolysolver {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/hla", mode: params.publishDirMode }
  
  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForPolysolver

  output:
    set idTumor, idNormal, target, file("${outputDir}/winners.hla.txt") into hlaOutput

  when: "polysolver" in tools && runSomatic
  
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
    ${outputDir}
  """
}

(bamsForConpair, bamFiles) = bamFiles.into(2)

// --- Run Conpair
process RunConpair {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/qc/conpair/${idTumor}_vs_${idNormal}", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForConpair
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set file("${idNormal}.pileup"), file("${idTumor}.pileup") into conpairPileup
    set file("${outPrefix}.concordance.txt"), file("${outPrefix}.contamination.txt") into conpairOutput

  when: !params.test

  script:
  outPrefix = "${idTumor}_vs_${idNormal}"
  gatkPath = "/usr/bin/GenomeAnalysisTK.jar"
  conpairPath = "/usr/bin/conpair"
  
  markersBed = ""
  markersTxt = ""
  if (params.genome == "GRCh37") {
    markersBed = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed"
    markersTxt = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
  }
  else {
    markersBed = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed"
    markersTxt = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
  }

  mem = 0
  if (params.mem_per_core) {
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
    --outpre=${outPrefix}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${idTumor}.pileup \
    --normal_pileup=${idNormal}.pileup \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${outPrefix}
  
  mv ${outPrefix}_concordance.txt ${outPrefix}.concordance.txt
  mv ${outPrefix}_contamination.txt ${outPrefix}.contamination.txt
  """
}

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

// create channel for MetaDataParser, using *purity.out from FACETS,

(facetsForLOHHLA, facetsForMetaDataParser) = facetsForLOHHLA.into(2)

(hlaOutputForLOHHLA, hlaOutput) = hlaOutput.into(2)

// *purity.out from FACETS, winners.hla.txt from POLYSOLVER, with the above

//apply *.groupTuple(by: [0,1,2]) in order to group the channel by idTumor, idNormal, and target

mergedChannelLOHHLA = bamsForLOHHLA.combine(hlaOutputForLOHHLA, by: [0,1,2]).combine(facetsForLOHHLA, by: [0,1,2])

// Run LOHHLA
process RunLOHHLA {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic/lohhla", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file("winners.hla.txt"), file("*_purity.out") from mergedChannelLOHHLA
    set file(hlaFasta), file(hlaDat) from Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat])

  output:
    set file("*HLAlossPrediction_CI.xls"), file("Figures/*.pdf") into lohhlaOutput

  when: tools.containsAll(["lohhla", "polysolver", "facets"]) && runSomatic

  script:
  """
  cat winners.hla.txt | tr "\t" "\n" | grep -v "HLA" > massaged.winners.hla.txt
  
  PURITY=\$(grep Purity *_purity.out | grep -oP "[0-9\\.]+")
  PLOIDY=\$(grep Ploidy *_purity.out | grep -oP "[0-9\\.]+")
  cat <(echo -e "tumorPurity\ttumorPloidy") <(echo -e "\$PURITY\t\$PLOIDY") > tumor_purity_ploidy.txt

  Rscript /lohhla/LOHHLAscript.R \
    --patientId ${idTumor}_vs_${idNormal} \
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

(mafFileForMafAnno, mafFileForMutSig, mafFile) = mafFile.into(3)

// --- Run Mutational Signatures, github.com/mskcc/mutation-signatures
process RunMutationSignatures {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll){ publishDir "${params.outDir}/somatic/mutation_signatures", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(maf) from mafFileForMutSig

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.mutsig.txt") into mutSigOutput

  when: tools.containsAll(["mutect2", "manta", "strelka2", "mutsig"]) && runSomatic

  script:
  outputPrefix = "${idTumor}_vs_${idNormal}"
  """
  python /mutation-signatures/main.py \
    /mutation-signatures/Stratton_signatures30.txt \
    ${outputPrefix}.somatic.maf \
    ${outputPrefix}.mutsig.txt
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
    
    return [idTumor, idNormal, target, purity_rdata, purity_cncf, hisens_cncf]
  }

//Formatting the channel to be grouped by idTumor, idNormal, and target

mafFileForMafAnno = mafFileForMafAnno.groupTuple(by: [0,1,2])

FacetsMafFileCombine = FacetsforMafAnno.combine(mafFileForMafAnno, by: [0,1,2])

// --- Do FACETS MAF annotation and post processing
process FacetsAnnotation {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(purity_rdata), file(purity_cncf), file(hisens_cncf), file(maf) from FacetsMafFileCombine

  output:
    set idTumor, idNormal, target, file("${outputPrefix}.facets.maf"), file("${outputPrefix}.armlevel.tsv"), file("${outputPrefix}.genelevel.tsv"), file("${outputPrefix}.genelevel_TSG_ManualReview.txt") into FacetsAnnotationOutputs

  when: tools.containsAll(["facets", "mutect2", "manta", "strelka2"]) && runSomatic

  script:
  mapFile = "${idTumor}_${idNormal}.map"
  outputPrefix = "${idTumor}_vs_${idNormal}"
  """
  echo "Tumor_Sample_Barcode\tRdata_filename" > ${mapFile}
  echo "${idTumor}\t${purity_rdata.fileName}" >> ${mapFile}

  /usr/bin/facets-suite/mafAnno.R \
    --facets_files ${mapFile} \
    --maf ${maf} \
    --out_maf ${outputPrefix}.facets.maf
  
  /usr/bin/facets-suite/geneLevel.R \
    --filenames ${hisens_cncf} \
    --outfile ${outputPrefix}.genelevel.tsv

  /usr/bin/facets-suite/armLevel.R \
    --filenames ${purity_cncf} \
    --outfile ${outputPrefix}.armlevel.tsv
  """
}

(mafFileForNeoantigen, facetsAnnotationForMetaData, FacetsAnnotationOutputs) = FacetsAnnotationOutputs.into(3)

//Formatting the channel to be: idTumor, idNormal, target, MAF

mafFileForNeoantigen = mafFileForNeoantigen.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def mafFile = item[3]
    def armLevel = item[4]
    def geneLevel = item[5]
    def tsg_manual_review = item[6]
    return [idTumor, idNormal, target, mafFile]
  }

hlaOutput = hlaOutput.combine(mafFileForNeoantigen, by: [0,1,2])

(hlaOutputForMetaDataParser, hlaOutput) = hlaOutput.into(2)

// --- Run neoantigen prediction pipeline
process RunNeoantigen {
  tag {idTumor + "_vs_" + idNormal}

  if (publishAll) { publishDir "${params.outDir}/somatic", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(polysolverFile), file(mafFile) from hlaOutput
    set file(neoantigenCDNA), file(neoantigenCDS) from Channel.value([
      referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS
    ])

  output:
    set idTumor, idNormal, target, file("${outputDir}/*") into neoantigenOut
    file("${idTumor}_vs_${idNormal}.all_neoantigen_predictions.txt") into NetMhcStatsOutput
    file("${outputDir}/*.maf") into NeoantigenMafOutput

  when: tools.containsAll(["neoantigen", "mutect2", "manta", "strelka2"]) && runSomatic

  script:
  outputDir = "neoantigen"
  tmpDir = "${outputDir}-tmp"
  tmpDirFullPath = "\$PWD/${tmpDir}/"  // must set full path to tmp directories for netMHC and netMHCpan to work; for some reason doesn't work with /scratch, so putting them in the process workspace
  """
  export TMPDIR=${tmpDirFullPath}
  mkdir -p ${tmpDir}
  chmod 777 ${tmpDir}

  python /usr/local/bin/neoantigen/neoantigen.py \
    --config_file /usr/local/bin/neoantigen/neoantigen-docker.config \
    --sample_id ${idTumor}_vs_${idNormal} \
    --hla_file ${polysolverFile} \
    --maf_file ${mafFile} \
    --output_dir ${outputDir}

  awk 'NR==1 {printf("%s\\t%s\\n", \$0, "sampleID")}  NR>1 {printf("%s\\t%s\\n", \$0, "${idTumor}_vs_${idNormal}") }' neoantigen/*.all_neoantigen_predictions.txt > ${idTumor}_vs_${idNormal}.all_neoantigen_predictions.txt

  """
}

// [idTumor, idNormal, target, armLevel]
facetsAnnotationForMetaData = facetsAnnotationForMetaData.map{
  item -> 
    def idTumor = item[0]
    def idNormal = item[1]
    def target = item[2]
    def mafFile = item[3]
    def armLevel = item[4]
    def geneLevel = item[5]
    def tsg_manual_review = item[6]
    return [idTumor, idNormal, target, armLevel]
  }

(mutsigMetaData, mutSigOutput) = mutSigOutput.into(2)

mergedChannelMetaDataParser = facetsForMetaDataParser.combine(facetsAnnotationForMetaData, by: [0,1,2]).combine(msiOutputForMetaData, by: [0,1,2]).combine(hlaOutputForMetaDataParser, by: [0,1,2]).combine(mutsigMetaData, by: [0,1,2]).unique()

// --- Generate sample-level metadata
process MetaDataParser {
  tag {idSample}

  if (publishAll) { publishDir "${params.outDir}/", mode: params.publishDirMode }
 
  input:
    set idTumor, idNormal, target, file(purityOut), file(armLevel), file(msifile), file(polysolverFile), file(mafFile), file(mutSigOutput) from mergedChannelMetaDataParser
    set file(idtCodingBed), file(agilentCodingBed), file(wgsCodingBed) from Channel.value([
      referenceMap.idtCodingBed, referenceMap.agilentCodingBed, referenceMap.wgsCodingBed
    ]) 

  output:
    file("*_metadata.tsv") into MetaDataOutputs

  when: runSomatic

  script:
  codingRegionsBed = ""
  if (target == "idt") {
    codingRegionsBed = "${idtCodingBed}"
  }
  else if (target == "agilent") {
    codingRegionsBed = "${agilentCodingBed}"
  }
  else if (target == "wgs") {
    codingRegionsBed = "${wgsCodingBed}"
  }
  """
  python3 /usr/bin/create_metadata_file.py \
    --sampleID ${idTumor}_vs_${idNormal} \
    --facetsPurity_out ${purityOut} \
    --facetsArmLevel  ${armLevel} \
    --MSIsensor_output ${msifile} \
    --mutational_signatures_output ${mutSigOutput} \
    --polysolver_output ${polysolverFile} \
    --MAF_input ${mafFile} \
    --coding_baits_BED ${codingRegionsBed}
  """
}

process SomaticAggregate {
 
  publishDir "${params.outDir}/somatic/", mode: params.publishDirMode

  input:
    file(netmhcCombinedFile) from NetMhcStatsOutput.collect()
    file(mafFile) from NeoantigenMafOutput.collect()
    file(mutsigFile) from mutSigOutput.collect()
    file(purityFiles) from FacetsPurity.collect()
    file(hisensFiles) from FacetsHisens.collect()
    file(purityHisensOutput) from FacetsPurityHisensOutput.collect()
    file(annotationFiles) from FacetsAnnotationOutputs.collect()
    file(dellyMantaVcf) from vcfDellyMantaMergedOutput.collect()
    file(metaDataFile) from MetaDataOutputs.collect()
    file(facetsOutputSubdirectories) from FacetsOutputSubdirectories.collect()

  output:
    file("merged.maf") into MafFileOutput
    file("merged_all_neoantigen_predictions.txt") into NetMhcChannel
    file("mutsig/*") into MutSigFilesOutput
    file("facets/*") into FacetsChannel
    set file("merged_hisens.cncf.txt"), file("merged_purity.cncf.txt"), file("merged_hisens.seg"), file("merged_purity.seg") into FacetsMergedChannel
    set file("merged_armlevel.tsv"), file("merged_armlevel.tsv"), file("merged_genelevel_TSG_ManualReview.txt"), file("merged_hisensPurity_out.txt") into FacetsAnnotationMergedChannel
    file("merged.vcf.gz") into VcfBedPeChannel
    file("merged_metadata.tsv") into MetaDataOutputChannel

  when: runSomatic
    
  script:
  """
  # Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  # Collect MAF files from neoantigen to maf_files/ and merge into one maf
  mkdir maf_files
  mv *.maf maf_files
  cat maf_files/*.maf | grep ^Hugo | head -n1 > merged.maf
  cat maf_files/*.maf | grep -Ev "^#|^Hugo" | sort -k5,5V -k6,6n >> merged.maf

  # Collect netmhc/netmhcpan combined files from neoantigen to netmhc_stats
  mkdir netmhc_stats
  mv *.all_neoantigen_predictions.txt netmhc_stats
  awk 'FNR==1 && NR!=1{next;}{print}' netmhc_stats/*.all_neoantigen_predictions.txt > merged_all_neoantigen_predictions.txt 

  # Collect mutsig output to mutsig/
  mkdir mutsig
  mv *.mutsig.txt mutsig/

  # Collect facets output to facets/
  mkdir facets
  mkdir facets/hisens
  mkdir facets/purity
  mkdir facets/hisensPurityOutput
  mv *purity.* facets/purity
  mv *hisens.* facets/hisens
  mv *_OUT.txt facets/hisensPurityOutput
  awk 'FNR==1 && NR!=1{next;}{print}' facets/hisens/*_hisens.cncf.txt > merged_hisens.cncf.txt
  awk 'FNR==1 && NR!=1{next;}{print}' facets/purity/*_purity.cncf.txt > merged_purity.cncf.txt
  awk 'FNR==1 && NR!=1{next;}{print}' facets/hisens/*_hisens.seg > merged_hisens.seg
  awk 'FNR==1 && NR!=1{next;}{print}' facets/purity/*_purity.seg > merged_purity.seg 
  awk 'FNR==1 && NR!=1{next;}{print}' facets/hisensPurityOutput/*_OUT.txt > merged_hisensPurity_out.txt  

  ## Move and merge FacetsAnnotation outputs
  mkdir facets/armLevel
  mkdir facets/geneLevel
  mkdir facets/manualReview
  mv *armlevel.tsv facets/armLevel
  mv *genelevel.tsv facets/geneLevel
  mv *genelevel_TSG_ManualReview.txt  facets/manualReview
  awk 'FNR==1 && NR!=1{next;}{print}' facets/armLevel/*armlevel.tsv > merged_armlevel.tsv
  awk 'FNR==1 && NR!=1{next;}{print}' facets/geneLevel/*genelevel.tsv > merged_genelevel.tsv
  awk 'FNR==1 && NR!=1{next;}{print}' facets/manualReview/*genelevel_TSG_ManualReview.txt > merged_genelevel_TSG_ManualReview.txt 

  ## Move all FACETS output subdirectories into /facets
  mv ${facetsOutputSubdirectories} facets/

  # Collect delly and manta vcf outputs into vcf_delly_manta/
  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done

  mkdir vcf_delly_manta
  mv *delly.manta.vcf.gz* vcf_delly_manta

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output merged.vcf.gz \
    vcf_delly_manta/*delly.manta.vcf.gz

  ## Collect metadata *tsv file into merged_metadata.tsv
  mkdir metadata
  mv *_metadata.tsv metadata 
  awk 'FNR==1 && NR!=1{next;}{print}' metadata/*_metadata.tsv > merged_metadata.tsv
  """
}

/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/

// GATK HaplotypeCaller
process GermlineRunHaplotypecaller {
  tag {idNormal + "@" + intervalBed.baseName}

  if (publishAll) { publishDir "${params.outDir}/germline/haplotypecaller", mode: params.publishDirMode }

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set id, target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set id, idTumor, idNormal, target, file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz"),
        file("${idNormal}_${intervalBed.baseName}.snps.filter.vcf.gz.tbi"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz"), file("${idNormal}_${intervalBed.baseName}.indels.filter.vcf.gz.tbi") into haplotypecallerOutput

  when: 'haplotypecaller' in tools && runGermline

  script:
  """
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

//Formatting the channel to be grouped by idTumor, idNormal, and target
// group by groupKey(key, intervalBed.size())
haplotypecallerOutput = haplotypecallerOutput.groupTuple()

// merge VCFs, GATK HaplotypeCaller
process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/haplotypecaller", mode: params.publishDirMode }

  input:
    set id, idTumor, idNormal, target, file(haplotypecallerSnpVcf), file(haplotypecallerSnpVcfIndex), file(haplotypecallerIndelVcf), file(haplotypecallerIndelVcfIndex) from haplotypecallerOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools && runGermline 


  script:
  outfile="${idNormal}.haplotypecaller.vcf.gz"
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
    --output ${idNormal}.haplotypecaller.vcf.gz

  tabix --preset vcf ${idNormal}.haplotypecaller.vcf.gz
  """
}

(bamsForMantaGermline, bamsForStrelkaGermline, bamFiles) = bamFiles.into(3)

// --- Run Manta, germline
process GermlineRunManta {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/manta", mode: params.publishDirMode }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForMantaGermline
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex
    ])
    set file(svCallingIncludeRegions), file(svCallingIncludeRegionsIndex) from Channel.value([
      referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex
    ])

  output:
    set idTumor, idNormal, target, file("Manta_${idNormal}.diploidSV.vcf.gz") into mantaOutputGermline mode flatten

  when: 'manta' in tools && runGermline

  // flag with --exome if exome
  script:
  options = ""
  if (assay == "wes") options = "--exome"
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

// --- Run Strelka2, germline
process GermlineRunStrelka2 {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/strelka2", mode: params.publishDirMode }

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForStrelkaGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])
    set file(idtTargets), file(agilentTargets), file(wgsIntervals),
    file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervalsIndex) from Channel.value([
      referenceMap.idtTargets, referenceMap.agilentTargets, referenceMap.wgsTargets,
      referenceMap.idtTargetsIndex, referenceMap.agilentTargetsIndex, referenceMap.wgsTargetsIndex
    ])
    
  output:
    set idTumor, idNormal, target, file("Strelka_${idNormal}_variants.vcf.gz"), file("Strelka_${idNormal}_variants.vcf.gz.tbi") into strelkaOutputGermline

  when: 'strelka2' in tools && runGermline
  
  script:
  options = ""
  intervals = wgsIntervals
  if (assay == "wes") {
    options = "--exome"
    if (target == 'agilent') intervals = agilentTargets
    if (target == 'idt') intervals = idtTargets
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

// Join HaploTypeCaller and Strelka outputs,  bcftools
haplotypecallerStrelkaChannel = haplotypecallerCombinedVcfOutput.combine(strelkaOutputGermline, by: [0,1,2])

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

mergedChannelVcfCombine = bamsForCombineChannel.combine(haplotypecallerStrelkaChannel, by: [0,1,2])

// --- Combine VCFs with germline calls from Haplotypecaller and Strelka2
process GermlineCombineChannel {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/haplotypecaller_strelka2_vcf", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, assay, file(bamTumor), file(baiTumor), file(haplotypecallercombinedVcf), file(haplotypecallercombinedVcfIndex), file(strelkaVcf), file(strelkaVcfIndex) from mergedChannelVcfCombine
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex,
    ])
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker, referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex
    ])
    set file(gnomadWesVcf), file(gnomadWesVcfIndex), file(gnomadWgsVcf), file(gnomadWgsVcfIndex) from Channel.value([
      referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,
      referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}.germline.vcf") into vcfMergedOutputGermline

  when: tools.containsAll(["strelka2", "haplotypecaller"]) && runGermline

  script:  
  isec_dir = "${idNormal}.isec"
  gnomad = gnomadWgsVcf
  if (target == 'wgs') {
    gnomad = gnomadWgsVcf
    gnomadCutoff = 'AF_popmax>0.02'
  }
  else {
    gnomad = gnomadWesVcf
    gnomadCutoff = 'non_cancer_AF_popmax>0.02'
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
    --exclude \"${gnomadCutoff}\" \
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

// Run vcf2maf on combined germline VCF, apply custom filters
process GermlineAnnotateMaf {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/mutations", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(vcfMerged) from vcfMergedOutputGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(vepCache), file(isoforms) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
      referenceMap.vepCache, referenceMap.isoforms
    ])

  output:
    file("${outputPrefix}.maf") into mafFileGermline

  when: tools.containsAll(["strelka2", "haplotypecaller"]) && runGermline

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


// --- Call germline SVs with Delly
svTypes = Channel.from("DUP", "BND", "DEL", "INS", "INV")
(bamsForDellyGermline, bamFiles) = bamFiles.into(2)

process GermlineDellyCall {
  tag {idNormal + '@' + svType}

  if (publishAll) { publishDir "${params.outDir}/germline/delly", mode: params.publishDirMode }

  input:
    each svType from svTypes
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamsForDellyGermline
    set file(genomeFile), file(genomeIndex), file(svCallingExcludeRegions) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}_${svType}.filter.bcf") into dellyFilterOutputGermline

  when: 'delly' in tools && runGermline

  script:
  """
  delly call \
    --svtype ${svType} \
    --genome ${genomeFile} \
    --exclude ${svCallingExcludeRegions} \
    --outfile ${idNormal}_${svType}.bcf \
    ${bamNormal}

  delly filter \
    --filter germline \
    --outfile ${idNormal}_${svType}.filter.bcf \
    ${idNormal}_${svType}.bcf
  """
}

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
// filter Delly & Manta via bcftools


// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)

// filter Delly & Manta via bcftools


dellyFilterOutputGermline = dellyFilterOutputGermline.groupTuple(by: [0,1,2], size: 5)

dellyMantaChannelGermline = dellyFilterOutputGermline.combine(mantaOutputGermline, by: [0,1,2])

// --- Merge Delly and Manta VCFs 
process GermlineMergeDellyAndManta {
  tag {idNormal}

  if (publishAll) { publishDir "${params.outDir}/germline/structural_variants", mode: params.publishDirMode }

  input:
    set idTumor, idNormal, target, file(dellyBcf), file(mantaVcf) from dellyMantaChannelGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}.delly.manta.vcf.gz"), file("${idNormal}.delly.manta.vcf.gz.tbi") into vcfFilterDellyMantaOutputGermline
    set file("${idNormal}.delly.manta.vcf.gz"), file("${idNormal}.delly.manta.vcf.gz.tbi") into germlineVcfBedPe

  when: tools.containsAll(["manta", "delly"]) && runGermline

  script:
  """ 
  for f in *.bcf
  do 
    bcftools view --output-type z \$f > \${f%.bcf}.delly.vcf.gz
  done

  for f in *.vcf.gz
  do
    tabix --preset vcf \$f
  done

  bcftools concat \
    --allow-overlaps \
    --output-type z \
    --output ${idNormal}.delly.manta.unfiltered.vcf.gz \
    *.delly.vcf.gz Manta_${idNormal}.diploidSV.vcf.gz

  tabix --preset vcf ${idNormal}.delly.manta.unfiltered.vcf.gz

  bcftools filter \
    --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.delly.manta.vcf.gz \
    ${idNormal}.delly.manta.unfiltered.vcf.gz 
    
  tabix --preset vcf ${idNormal}.delly.manta.vcf.gz
  """
}

// --- Aggregate per-sample germline data
process GermlineAggregate {
 
  publishDir "${params.outDir}/germline/", mode: params.publishDirMode

  input:
    file(mafFile) from mafFileGermline.collect()
    file(dellyMantaVcf) from germlineVcfBedPe.collect()

  output:
    file("merged.maf") into GermlineMafFileOutput
    file("merged.vcf.gz") into GermlineVcfBedPeChannel
  
  when: runGermline

  script:
  """
  # Making a temp directory that is needed for some reason...
  mkdir tmp
  TMPDIR=./tmp

  # Collect MAF files from neoantigen to maf_files/ and merge into one maf
  mkdir maf_files
  mv *.maf maf_files
  cat maf_files/*.maf | grep ^Hugo | head -n1 > merged.maf
  cat maf_files/*.maf | grep -Ev "^#|^Hugo" | sort -k5,5V -k6,6n >> merged.maf

  # Collect delly and manta vcf outputs into vcf_delly_manta/
  mkdir vcf_delly_manta
  mv  *.delly.manta.vcf.gz* vcf_delly_manta

  bcftools merge \
    --force-samples \
    --merge none \
    --output-type z \
    --output merged.vcf.gz \
    vcf_delly_manta/*.delly.manta.vcf.gz
  """
}

/*
================================================================================
=                              AUXILLIARY FUNCTIONS                            =
================================================================================
*/

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  result_array = [
    'dbsnp' : checkParamReturnFile("dbsnp"),
    'dbsnpIndex' : checkParamReturnFile("dbsnpIndex"),
    // genome reference dictionary
    'genomeDict' : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile' : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex' : checkParamReturnFile("genomeIndex"),
    // BWA index files
    'bwaIndex' : checkParamReturnFile("bwaIndex"), 
    // VCFs with known indels (such as 1000 Genomes, Mill’s gold standard)
    'knownIndels' : checkParamReturnFile("knownIndels"),
    'knownIndelsIndex' : checkParamReturnFile("knownIndelsIndex"),
    'msiSensorList' : checkParamReturnFile("msiSensorList"),
    'svCallingExcludeRegions' : checkParamReturnFile("svCallingExcludeRegions"),
    'svCallingIncludeRegions' : checkParamReturnFile("svCallingIncludeRegions"),
    'svCallingIncludeRegionsIndex' : checkParamReturnFile("svCallingIncludeRegionsIndex"),
    // Target and Bait BED files
    'idtTargets' : checkParamReturnFile("idtTargets"),
    'idtTargetsIndex' : checkParamReturnFile("idtTargetsIndex"),
    'idtTargetsList' : checkParamReturnFile("idtTargetsList"),  
    'idtBaitsList' : checkParamReturnFile("idtBaitsList"), 
    'agilentTargets' : checkParamReturnFile("agilentTargets"),
    'agilentTargetsIndex' : checkParamReturnFile("agilentTargetsIndex"),
    'agilentTargetsList' : checkParamReturnFile("agilentTargetsList"),  
    'agilentBaitsList' : checkParamReturnFile("agilentBaitsList"), 
    'wgsTargets' : checkParamReturnFile("wgsTargets"),
    'wgsTargetsIndex' : checkParamReturnFile("wgsTargetsIndex")
  ]

  if (!params.test) {
    result_array << ['vepCache' : checkParamReturnFile("vepCache")]
    // for SNP Pileup
    result_array << ['facetsVcf' : checkParamReturnFile("facetsVcf")]
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
    result_array << ['agilentCodingBed' : checkParamReturnFile("agilentCodingBed")]    
    result_array << ['wgsCodingBed' : checkParamReturnFile("wgsCodingBed")]  
  }
  return result_array
}

def debug(channel) {
  channel.subscribe { Object obj ->
    println "DEBUG: ${obj.toString()};"
  }
}

def extractPairing(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    [row.TUMOR_ID, row.NORMAL_ID]
  }
}

def extractFastq(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    checkNumberOfItem(row, 6)
    def idSample = row.SAMPLE
    def lane = row.LANE
    def assayValue = row.ASSAY
    def targetFile = row.TARGET
    def fastqFile1 = returnFile(row.FASTQ_PE1)
    def sizeFastqFile1 = fastqFile1.size()
    def fastqFile2 = returnFile(row.FASTQ_PE2)
    def sizeFastqFile2 = fastqFile2.size()

    def assay = assayValue.toLowerCase() //standardize genome/wgs/WGS to wgs, exome/wes/WES to wes

    if ((assay == "genome") || (assay == "wgs")) {
      assay = "wgs"
    }
    if ((assay == "exome") || (assay == "wes")) {
      assay = "wes"
    }

    checkFileExtension(fastqFile1,".fastq.gz")
    checkFileExtension(fastqFile2,".fastq.gz")

    [idSample, lane, fastqFile1, sizeFastqFile1, fastqFile2, sizeFastqFile2, assay, targetFile]
  }
}

def extractBAM(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    checkNumberOfItem(row, 6)
    def idTumor = row.TUMOR_ID
    def idNormal = row.NORMAL_ID
    def assay = row.ASSAY
    def target = row.TARGET
    def bamTumor = returnFile(row.TUMOR_BAM)
    // check if using bamTumor.bai or bamTumor.bam.bai
    def baiTumor = returnFile(validateBamIndexFormat(row.TUMOR_BAM))
    // def sizeTumorBamFile = tumorBamFile.size()
    def bamNormal = returnFile(row.NORMAL_BAM)
    def baiNormal = returnFile(validateBamIndexFormat(row.NORMAL_BAM))
    // def sizeNormalBamFile = normalBamFile.size()

    [assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)]
  }
}

// Check which format of BAM index used, input 'it' as BAM file 'bamTumor.bam'
def validateBamIndexFormat(it) {
  bamFilename = it.take(it.lastIndexOf('.'))
  // Check BAM index extension
  if (file(bamFilename + ".bai").exists()){
    return(file("${bamFilename}.bai"))
  } else if (file(bamFilename + ".bam.bai").exists()){
    return(file("${bamFilename}.bam.bai"))
  } else {
    println "ERROR: Cannot find BAM indices for ${it}. Please index BAMs in the same directory with 'samtools index' and re-run the pipeline."
    exit 1
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

def check_for_duplicated_rows(pairingFilePath) {
  def entries = []
  file( pairingFilePath ).eachLine { line ->
    if (!line.isEmpty()){
      entries << line
    }
  }
  return entries.toSet().size() == entries.size()
}

def check_for_mixed_assay(mappingFilePath) {
  def wgs = false
  def wes = false
  file( mappingFilePath ).eachLine { line ->
    currentLine = line.toLowerCase()
    if (currentLine.contains('\tgenome\t') || currentLine.contains('\twgs\t')) {
      wgs = true
    }
    if (currentLine.contains('\texome\t') || currentLine.contains('\twes\t')) {
      wes = true
    }
  return !(wgs && wes)
  }
}

// check lane names are unique in input mapping *tsv 
def checkForUniqueSampleLanes(inputFilename) {
  def totalList = []
  // parse tsv
  file(inputFilename).eachLine { line ->
      if (!line.isEmpty()){
          def (sample, lane, assay, target, fastqpe1, fastqpe2) = line.split(/\t/)
          totalList << sample + "_" + lane
      }
  }
  // remove header 'SAMPLE_LANE'
  totalList.removeAll{ it == 'SAMPLE_LANE'} 
  return totalList.size() == totalList.unique().size()
}
