#!/usr/bin/env nextflow

/*
================================================================================
--------------------------------------------------------------------------------
 Processes overview
 - AlignReads - Map reads with BWA mem output SAM
 - SortBAM - Sort BAM with samtools
 - MergeBam - Merge BAM for the same samples from different lanes
 - MarkDuplicates - Mark Duplicates with GATK4
 - CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
 - RecalibrateBam - Recalibrate Bam with PrintReads
*/

/*
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

if (params.mapping) mappingPath = params.mapping
if (params.pairing) pairingPath = params.pairing

outname = params.outname

referenceMap = defineReferenceMap()

fastqFiles = Channel.empty()

mappingFile = file(mappingPath)
pairingfile = file(pairingPath)

pairingT = extractPairing(pairingfile)

fastqFiles = extractFastq(mappingFile)

runGermline = params.germline
runSomatic = params.somatic

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/


fastqFiles.groupTuple(by:[0]).map { key, lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets -> tuple( groupKey(key, lanes.size()), lanes, files_pe1, files_pe1_size, files_pe2, files_pe2_size, assays, targets) }.set { groupedFastqs }

groupedFastqs.into { groupedFastqsDebug; fastPFiles; fastqFiles }
fastPFiles = fastPFiles.transpose()
fastqFiles = fastqFiles.transpose()

// AlignReads - Map reads with BWA mem output SAM

process AlignReads {
  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  input:
    set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assay, targetFile from fastqFiles
    set file(genomeFile), file(bwaIndex) from Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex])

  output:
    set idSample, lane, file("${lane}.bam"), assay, targetFile into (unsortedBam)

  script:
    readGroup = "@RG\\tID:${lane}\\tSM:${idSample}\\tLB:${idSample}\\tPL:Illumina"
    
  """
  set -e
  set -o pipefail
  bwa mem -R \"${readGroup}\" -t ${task.cpus} -M ${genomeFile} ${fastqFile1} ${fastqFile2} | samtools view -Sb - > ${lane}.bam
  """
}

// SortBAM - Sort unsorted BAM with samtools, 'samtools sort'

process SortBAM {
  tag {lane}

  input:
    set idSample, lane, file("${lane}.bam"), assay, targetFile from unsortedBam

  output:
    set idSample, lane, file("${lane}.sorted.bam"), assay, targetFile into (sortedBam, sortedBamDebug)

  script:
  // Refactor when https://github.com/nextflow-io/nextflow/pull/1035 is merged
  if(params.mem_per_core) { 
    mem = task.memory.toString().split(" ")[0].toInteger() - 1 
  }
  else {
    mem = (task.memory.toString().split(" ")[0].toInteger()/task.cpus).toInteger() - 1
  }
  """
  samtools sort -m ${mem}G -@ ${task.cpus} -o ${lane}.sorted.bam ${lane}.bam
  """
}

sortedBam.groupTuple().set { groupedBam }
groupedBam.into { groupedBamDebug; groupedBam }

// MergeBams

process MergeBams {
  tag {idSample}

  input:
    set idSample, lane, file(bam), assay, targetFile from groupedBam

  output:
    set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile into (mergedBam, mergedBamDebug)

  script:
  """
  samtools merge --threads ${task.cpus} ${idSample}.merged.bam ${bam.join(" ")}
  """
}

// GATK MarkDuplicates

process MarkDuplicates {
  tag {idSample}

   publishDir "${params.outDir}/MarkDup/${idSample}", mode: params.publishDirMode

  input:
    set idSample, lane, file("${idSample}.merged.bam"), assay, targetFile from mergedBam

  output:
    set file("${idSample}.md.bam"), file("${idSample}.md.bai"), idSample, lane, assay, targetFile into duplicateMarkedBams
    set idSample, val("${idSample}.md.bam"), val("${idSample}.md.bai"), assay, targetFile into markDuplicatesTSV
    file ("${idSample}.bam.metrics") into markDuplicatesReport

  script:
  """
  gatk MarkDuplicates --java-options ${params.markdup_java_options}  \
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
    bam, bai, idSample, lane, assay, targetFile ->
    tag = bam.baseName.tokenize('.')[0]
    [idSample, bam, bai, assay, targetFile]
}

(mdBam, mdBamToJoin, mdDebug) = duplicateMarkedBams.into(3)

// GATK BaseRecalibrator , CreateRecalibrationTable

process CreateRecalibrationTable {
  tag {idSample}

  input:
    set idSample, file(bam), file(bai), assay, targetFile from mdBam 

    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(knownIndels), file(knownIndelsIndex)  from Channel.value([
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
  known = knownIndels.collect{ "--known-sites ${it}" }.join(' ')

  """
  gatk BaseRecalibrator \
    --tmp-dir /tmp \
    --reference ${genomeFile} \
    --known-sites ${dbsnp} \
    ${known} \
    --verbosity INFO \
    --input ${bam} \
    --output ${idSample}.recal.table
  """
}

recalibrationTable = mdBamToJoin.join(recalibrationTable, by:[0])

(recalibrationTable, recalibrationTableDebug) = recalibrationTable.into(2)

// GATK ApplyBQSR, RecalibrateBAM

process RecalibrateBam {
  tag {idSample}

  publishDir "${params.outDir}/BQSR/${idSample}", mode: params.publishDirMode

  input:
    set idSample, file(bam), file(bai), assay, targetFile, file(recalibrationReport) from recalibrationTable

    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict 
    ])

  output:
    set idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bai"), assay, targetFile into recalibratedBam, recalibratedBamForStats, recalibratedBamForOutput, recalibratedBamForOutput2, recalibratedBamForDebug
    set idSample, val("${idSample}.recal.bam"), val("${idSample}.recal.bai"), assay, targetFile into recalibratedBamTSV
    set idSample into currentSample
    set file("${idSample}.recal.bam") into currentBam
    set file("${idSample}.recal.bai") into currentBai
    set assay into assays
    set targetFile into targets

  script:
  """
  gatk ApplyBQSR \
    --reference ${genomeFile} \
    --create-output-bam-index true \
    --bqsr-recal-file ${recalibrationReport} \
    --input ${bam} \
    --output ${idSample}.recal.bam
  """
}

recalibratedBamForOutput.combine(pairingT)
                        .filter { item -> // only keep combinations where sample is same as tumor pair sample
                          def sampleID = item[0]
                          def sampleBam = item[1]
                          def sampleBai = item[2]
                          def assay = item[3]
                          def target = item[4]
                          def tumorID = item[5]
                          def normalID = item[6]
                          sampleID == tumorID
                        }.map { item -> // re-order the elements
                          def sampleID = item[0]
                          def sampleBam = item[1]
                          def sampleBai = item[2]
                          def assay = item[3][0]
                          def target = item[4][0]
                          def tumorID = item[5]
                          def normalID = item[6]
                          def tumorBam = sampleBam
                          def tumorBai = sampleBai

                          return [ assay, target, tumorID, normalID, tumorBam, tumorBai ]
                        }.combine(recalibratedBamForOutput2)
                        .filter { item ->
                          def assay = item[0]
                          def target = item[1]
                          def tumorID = item[2]
                          def normalID = item[3]
                          def tumorBam = item[4]
                          def tumorBai = item[5]
                          def sampleID = item[6]
                          def normalBam = item[7]
                          def normalBai = item[8]
                          sampleID == normalID
                        }.map { item -> // re-order the elements
                          def assay = item[0]
                          def target = item[1]
                          def tumorID = item[2]
                          def normalID = item[3]
                          def tumorBam = item[4]
                          def tumorBai = item[5]
                          def sampleID = item[6]
                          def normalBam = item[7]
                          def normalBai = item[8]

                          return [ assay, target, tumorID, normalID, tumorBam, normalBam, tumorBai, normalBai ]
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

// FastP - FastP on lane pairs, R1/R2

process FastP {

  tag {lane}   // The tag directive allows you to associate each process executions with a custom label

  publishDir "${params.outDir}/FastP/${idSample}", mode: params.publishDirMode

  input:
    set idSample, lane, file(fastqFile1), sizeFastqFile1, file(fastqFile2), sizeFastqFile2, assays, targetFiles from fastPFiles
    
  output:
    file("*.html") into fastPResults 

  script:
  """
  fastp -h ${lane}.html -i ${fastqFile1} -I ${fastqFile2}
  """
}


ignore_read_groups = Channel.from( true , false )

// Alfred, BAM QC

process Alfred {
  tag {idSample}

  publishDir "${params.outDir}/Alfred/${idSample}", mode: params.publishDirMode
  
  input:
    each ignore_rg from ignore_read_groups
    set idSample, file(bam), file(bai), assay, targetFile from recalibratedBam

    file(genomeFile) from Channel.value([
      referenceMap.genomeFile
    ])

  output:
    set ignore_rg, idSample, file("*.tsv.gz"), file("*.tsv.gz.pdf") into bamsQCStats

  script:
  def ignore = ignore_rg ? "--ignore" : ''
  def outfile = ignore_rg ? "${idSample}.alfred.tsv.gz" : "${idSample}.alfred.RG.tsv.gz"
  """
  alfred qc --reference ${genomeFile} ${ignore} --outfile ${outfile} ${bam} && Rscript /opt/alfred/scripts/stats.R ${outfile}
  """

}

(sampleIdsForIntervalBeds, bamFiles) = bamFiles.into(2)

// GATK SplitIntervals, CreateScatteredIntervals

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

( bamsForIntervals, bamFiles ) = bamFiles.into(2)

//Associating interval_list files with BAM files, putting them into one channel
agilentIList = agilentIntervals.map{ n -> [ n, "agilent" ] }
idtIList = idtIntervals.map{ n -> [ n, "idt" ] }
wgsIList = wgsIntervals.map{ n -> [ n, "wgs" ] }

( aBamList, iBamList, wBamList ) = bamsForIntervals.into(3)

aMergedChannel = aBamList.combine(agilentIList, by: 1).unique() 
bMergedChannel = iBamList.combine(idtIList, by: 1).unique() 
wMergedChannel = wBamList.combine(wgsIList, by: 1).unique() 

// These will go into mutect2 and haplotypecaller
( mergedChannelSomatic, mergedChannelGermline  ) = aMergedChannel.concat( bMergedChannel, wMergedChannel).into(2) // { mergedChannelSomatic, mergedChannelGermline }

/*
================================================================================
=                                SOMATIC PIPELINE                              =
================================================================================
*/

// parse --tools parameter for downstream 'when' conditionals, e.g. when: `` 'delly ' in tools
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

// --- Run Mutect2

process RunMutect2 {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2", mode: params.publishDirMode

  input:
    // Order has to be target, assay, etc. because the channel gets rearranged on ".combine"
    set target, assay, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal), file(intervalBed) from mergedChannelSomatic 
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz"), file("${idTumor}_vs_${idNormal}_${intervalBed.baseName}.vcf.gz.tbi") into mutect2Output

  when: 'mutect2' in tools && runSomatic

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

// GATK FilterMutectCalls, RunMutect2Filter

process RunMutect2Filter {
  tag {idTumor + "_vs_" + idNormal + '_' + mutect2Vcf.baseName}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/mutect2_filter", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(mutect2Vcf), file(mutect2VcfIndex) from mutect2Output

  output:
    set idTumor, idNormal, target, file("*filtered.vcf.gz"), file("*filtered.vcf.gz.tbi"), file("*Mutect2FilteringStats.tsv") into forMutect2Combine mode flatten

  when: 'mutect2' in tools && runSomatic

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

// Combine Mutect2 VCFs, bcftools

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

  when: 'mutect2' in tools && runSomatic

  script:
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

  when: 'manta' in tools && runSomatic

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

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)

mantaOutput = mantaOutput.groupTuple(by: [0,1,2])

dellyFilterOutput = dellyFilterOutput.groupTuple(by: [0,1,2])

dellyMantaCombineChannel = dellyFilterOutput.combine(mantaOutput, by: [0,1,2]).unique()

// --- Process Delly and Manta VCFs 

(sampleIdsForDellyMantaMerge, bamFiles) = bamFiles.into(2)

// Merge VCFs, Delly and Manta

process SomaticMergeDellyAndManta {
  tag {idTumor + "_vs_" + idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/vcf_merged_output", mode: params.publishDirMode

  input:
    set idTumor, idNormal, target, file(dellyBcfs), file(mantaFile) from dellyMantaCombineChannel

  output:
    file("*filtered.merge.vcf") into vcfDellyMantaMergedOutput

  when: 'manta' in tools && 'delly' in tools && runSomatic

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

  when: 'manta' in tools && 'strelka2' in tools && runSomatic

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

  when: 'manta' in tools && 'strelka2' in tools && runSomatic

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

// Combined Somatic VCFs

process SomaticCombineChannel {
  tag {idTumor + "_vs_" + idNormal}

  input:
    set idTumor, idNormal, target, file(mutect2combinedVCF), file(mutect2combinedVCFIndex), file(strelkaVCF), file(strelkaVCFIndex) from mutectStrelkaChannel
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

  output:
    set idTumor, idNormal, target, file("${idTumor}.union.pass.vcf") into vcfMergedOutput

  when: 'manta' in tools && 'strelka2' in tools && 'mutect2' in tools && runSomatic

  script:
  isec_dir = "${idTumor}.isec"
  pon = wgsPoN
  if (target != 'wgs') {
    pon = exomePoN
  }
  """
  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">\n##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">\n##INFO=<ID=Strelka2Fail,Number=0,Type=Flag,Description=\"Variant was called failed by Strelka2\">" > vcf.header
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
    --mark-sites \"+Strelka2FAIL\" \
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
    --columns +FILTER,+FORMAT,Strelka2FAIL \
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
    --header-lines vcf.pon.header \
    --annotations ${pon} \
    --columns PoN:=AC \
    --output-type z \
    --output ${idTumor}.union.pon.vcf.gz \
    ${idTumor}.union.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idTumor}.union.pass.vcf \
    ${idTumor}.union.pon.vcf.gz
  """
}

// run VCF2MAF, somatic

process SomaticRunVcf2Maf {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${ params.outDir }/${idTumor}_vs_${idNormal}/somatic_variants/mutations", mode: params.publishDirMode

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
    file("*.maf") into mafFile

  // when: "mutect2" in tools && "manta" in tools && "strelka2" in tools && runSomatic 

  script:
  outfile="${vcfMerged}".replaceFirst(".vcf", ".maf")
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
    --retain-info MuTect2,Strelka2,RepeatMasker,EncodeDacMapability,PoN \
    --custom-enst ${isoforms} \
    --output-maf ${outfile} \
    --filter-vcf 0
  """
}

(bamFilesForMsiSensor, bamFiles) = bamFiles.into(2)

// MSI Sensor

process RunMsiSensor {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/msisensor", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal) from bamFilesForMsiSensor
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(msiSensorList) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.msiSensorList
    ])

  output:
    file("${output_prefix}*") into msiOutput 

  when: "msisensor" in tools && runSomatic

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

  when: 'facets' in tools && runSomatic

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

// FACETS

process DoFacets {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/facets", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(snpPileupFile) from SnpPileup

  output:
    file("*.*") into FacetsOutput

  when: 'facets' in tools && runSomatic

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

(bamsForHlaPolysolver, bamFiles) = bamFiles.into(2)

// Run HLA Polysolver

process RunHlaPolysolver {
  tag {idTumor + "_vs_" + idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/somatic_variants/hla_polysolver", mode: params.publishDirMode

  input:
    set assay, target, idTumor, idNormal, file(bamTumor), file(bamNormal), file(baiTumor), file(baiNormal)  from bamsForHlaPolysolver

  output:
    file("output/*") into hlaOutput

  when: "hla" in tools && runSomatic
  
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

  when: 'conpair' in tools && runSomatic

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

/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/

// GATK HaplotypeCaller

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
    set idTumor, idNormal, target, file("*.vcf.gz"), file("*.vcf.gz.tbi") into haplotypecallerOutput mode flatten

  when: 'haplotypecaller' in tools && runGermline

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

//Formatting the channel to be keyed by idTumor, idNormal, and target
haplotypecallerOutput = haplotypecallerOutput.groupTuple(by: [0,1,2])

// merge VCFs, GATK HaplotypeCaller

process GermlineCombineHaplotypecallerVcf {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/haplotypecaller"

  input:
    set idTumor, idNormal, target, file(haplotypecallerVcf), file(haplotypecallerVcfIndex) from haplotypecallerOutput
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("${outfile}"), file("${outfile}.tbi") into haplotypecallerCombinedVcfOutput

  when: 'haplotypecaller' in tools && runGermline 

  script:
  outfile="${idNormal}.haplotypecaller.filtered.vcf.gz"
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
    set idTumor, idNormal, target, file("*.vcf.gz") into mantaOutputGermline mode flatten

  when: 'manta' in tools && runGermline

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
    set file(idtTargetsIndex), file(agilentTargetsIndex), file(wgsIntervals) from Channel.value([
      referenceMap.idtTargetsIndex,
      referenceMap.agilentTargetsIndex,
      referenceMap.wgsTargetsIndex
      ])

  output:
    set idTumor, idNormal, target, file("Strelka_${idNormal}_variants.vcf.gz"), file("Strelka_${idNormal}_variants.vcf.gz.tbi") into strelkaOutputGermline

  when: 'strelka2' in tools && runGermline
  
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

// Join HaploTypeCaller and Strelka outputs,  bcftools

hcv = haplotypecallerCombinedVcfOutput.groupTuple(by: [0,1,2])

haplotypecallerStrelkaChannel = hcv.combine(strelkaOutputGermline, by: [0,1,2]).unique()

process GermlineCombineChannel {
  tag {idNormal}

  input:
    set idTumor, idNormal, target, file(haplotypecallercombinedVCF), file(haplotypecallercombinedVCFIndex), file(strelkaVCF), file(strelkaVCFIndex) from haplotypecallerStrelkaChannel
    set file(repeatMasker), file(repeatMaskerIndex), file(mapabilityBlacklist), file(mapabilityBlacklistIndex) from Channel.value([
      referenceMap.repeatMasker,
      referenceMap.repeatMaskerIndex,
      referenceMap.mapabilityBlacklist,
      referenceMap.mapabilityBlacklistIndex
    ])

  output:
    set idTumor, idNormal, target, file("${idNormal}.union.pass.vcf") into vcfMergedOutputGermline

  when: 'strelka2' in tools && 'haplotypecaller' in tools && runGermline

  script:  
  isec_dir = "${idNormal}.isec"
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
    --output ${idNormal}.union.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${idNormal}.union.pass.vcf \
    ${idNormal}.union.vcf
  """
}

// vcf2maf, germline calls

process GermlineRunVcf2Maf {
  tag {idNormal}

  publishDir "${ params.outDir }/${idTumor}_vs_${idNormal}/germline_variants/mutations"

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
    set idTumor, idNormal, target, file("*.maf") into mafFileGermline

  when: "strelka2" in tools && "haplotypecaller" in tools && runGermline

  // both tumor-id and normal-id flags are set to idNormal since we're not processing the tumor in germline.nf
  script:
  outfile="${vcfMerged}".replaceFirst(".vcf", ".maf")
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

  when: 'delly' in tools && runGermline

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

// filter germline Delly calls, bcftools

process GermlineDellyFilter {
  tag {idNormal + '_' + svType}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/delly"

  input:
    set idTumor, idNormal, target, svType, file(dellyBcf), file(dellyBcfIndex) from dellyCallOutputGermline


  output:
    set idTumor, idNormal, target, file("*.filter.bcf") into dellyFilterOutputGermline


  when: 'delly' in tools && runGermline

  outfile="${dellyBcf}".replaceFirst(".bcf",".filter.bcf")

  script:
  """
  delly filter \
    --filter germline \
    --outfile ${outfile} \
    ${dellyBcf}
  """
}

// Put manta output and delly output into the same channel so they can be processed together in the group key
// that they came in with i.e. (`idTumor`, `idNormal`, and `target`)

mantaOutputGermline = mantaOutputGermline.groupTuple(by: [0,1,2])
dellyFilterOutputGermline = dellyFilterOutputGermline.groupTuple(by: [0,1,2])

dellyMantaChannelGermline = dellyFilterOutputGermline.combine(mantaOutputGermline, by: [0,1,2]).unique()

process GermlineMergeDellyAndManta {
  tag {idNormal}

  //publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/vcf_merged_output"

  input:
    set idTumor, idNormal, target, file(dellyBcf), file(mantaVcf) from dellyMantaChannelGermline

  output:
    set idTumor, idNormal, target, file("*.merge.vcf.gz") into vcfDellyMantaMergedOutputGermline

  when: 'manta' in tools && 'delly' in tools && runGermline

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

// filter Delly & Manta via bcftools

process GermlineRunBcfToolsFilterOnDellyManta {
  tag {idNormal}

  publishDir "${params.outDir}/${idTumor}_vs_${idNormal}/germline_variants/mutations"

  input:
    set idTumor, idNormal, target, file(vcf) from vcfDellyMantaMergedOutputGermline
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idTumor, idNormal, target, file("*filtered.vcf.gz") into vcfFilterDellyMantaOutputGermline

  when: "manta" in tools && "delly" in tools && runGermline

  script:
  outfile = "${vcf}".replaceFirst('vcf.gz', 'filtered.vcf.gz')
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
  result_array =  [
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

def debug(channel) {
  channel.subscribe { Object obj ->
    println "DEBUG: ${obj.toString()};"
  }
}

def extractPairing(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    [row.NORMAL_ID, row.TUMOR_ID]
  }
}

def extractFastq(tsvFile) {
  Channel.from(tsvFile)
  .splitCsv(sep: '\t', header: true)
  .map { row ->
    checkNumberOfItem(row, 6)
    def idSample = row.SAMPLE
    def lane = row.LANE
    def assay = row.ASSAY
    def targetFile = row.TARGET
    def fastqFile1 = returnFile(row.FASTQ_PE1)
    def sizeFastqFile1 = fastqFile1.size()
    def fastqFile2 = returnFile(row.FASTQ_PE2)
    def sizeFastqFile2 = fastqFile2.size()

    checkFileExtension(fastqFile1,".fastq.gz")
    checkFileExtension(fastqFile2,".fastq.gz")

    [idSample, lane, fastqFile1, sizeFastqFile1, fastqFile2, sizeFastqFile2, assay, targetFile]
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
