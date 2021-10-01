#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println 'ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)'
  exit 1
}

// User-set runtime parameters
outDir     = file(params.outDir).toAbsolutePath()
outname = params.outname
runGermline = params.germline
runSomatic = params.somatic
runQC = params.QC
runAggregate = params.aggregate
runConpairAll = false
wallTimeExitCode = params.wallTimeExitCode ? params.wallTimeExitCode.split(',').collect { it.trim().toLowerCase() } : []
multiqcWesConfig = workflow.projectDir + '/lib/multiqc_config/exome_multiqc_config.yaml'
multiqcWgsConfig = workflow.projectDir + '/lib/multiqc_config/wgs_multiqc_config.yaml'
multiqcTempoLogo = workflow.projectDir + '/docs/tempoLogo.png'
epochMap = [:]
startEpoch = new Date().getTime()
limitInputLines = 0
chunkSizeLimit = params.chunkSizeLimit

include { defineReferenceMap; loadTargetReferences } from './modules/local/define_maps'
include { touchInputs; watchMapping; watchBamMapping; watchPairing; watchAggregateWithPath; watchAggregate } from './modules/local/watch_inputs'

if (params.watch == true) {
  touchInputs(chunkSizeLimit, startEpoch, epochMap)
}

//Alignment
include { SplitLanesR1; SplitLanesR2 } from './modules/process/Alignment/SplitLanes' 
include { AlignReads }                 from './modules/process/Alignment/AlignReads'
include { MergeBamsAndMarkDuplicates } from './modules/process/Alignment/MergeBamsAndMarkDuplicates'
include { RunBQSR }                    from './modules/process/Alignment/RunBQSR' 
include { CrossValidateSamples }       from './modules/process/QC/SampleValidation'

//Somatic
include { CreateScatteredIntervals }   from './modules/process/Somatic/CreateScatteredIntervals'
include { RunMutect2 }                 from './modules/process/Somatic/RunMutect2'
include { SomaticCombineMutect2Vcf }   from './modules/process/Somatic/SomaticCombineMutect2Vcf' 
include { SomaticDellyCall }           from './modules/process/Somatic/SomaticDellyCall' 
include { SomaticRunManta }            from './modules/process/Somatic/SomaticRunManta' 
include { SomaticMergeDellyAndManta }  from './modules/process/Somatic/SomaticMergeDellyAndManta' 
include { SomaticRunStrelka2 }         from './modules/process/Somatic/SomaticRunStrelka2' 
include { SomaticCombineChannel }      from './modules/process/Somatic/SomaticCombineChannel' 
include { SomaticAnnotateMaf }         from './modules/process/Somatic/SomaticAnnotateMaf' 
include { RunMutationSignatures }      from './modules/process/Somatic/RunMutationSignatures'
include { DoFacets }                   from './modules/process/Somatic/DoFacets' 
include { DoFacetsPreviewQC }          from './modules/process/Somatic/DoFacetsPreviewQC' 
include { RunPolysolver }              from './modules/process/Somatic/RunPolysolver'
include { RunLOHHLA }                  from './modules/process/Somatic/RunLOHHLA' 
include { RunNeoantigen }              from './modules/process/Somatic/RunNeoantigen' 
include { SomaticFacetsAnnotation }    from './modules/process/Somatic/SomaticFacetsAnnotation' 
include { RunMsiSensor }               from './modules/process/Somatic/RunMsiSensor'
include { MetaDataParser }             from './modules/process/Somatic/MetaDataParser' 

//Germline
include { GermlineRunHaplotypecaller }        from './modules/process/Germline/GermlineRunHaplotypecaller'
include { GermlineCombineHaplotypecallerVcf } from './modules/process/Germline/GermlineCombineHaplotypecallerVcf' 
include { GermlineRunStrelka2 }               from './modules/process/Germline/GermlineRunStrelka2' 
include { GermlineCombineChannel }            from './modules/process/Germline/GermlineCombineChannel' 
include { GermlineAnnotateMaf }               from './modules/process/Germline/GermlineAnnotateMaf' 
include { GermlineFacetsAnnotation }          from './modules/process/Germline/GermlineFacetsAnnotation' 
include { GermlineDellyCall }                 from './modules/process/Germline/GermlineDellyCall' 
include { GermlineRunManta }                  from './modules/process/Germline/GermlineRunManta' 
include { GermlineMergeDellyAndManta }        from './modules/process/Germline/GermlineMergeDellyAndManta'

//QC
include { QcCollectHsMetrics }                 from './modules/process/QC/QcCollectHsMetrics' 
include { QcQualimap }                         from './modules/process/QC/QcQualimap' 
include { QcAlfred }                           from './modules/process/QC/QcAlfred'
include { SampleRunMultiQC }                   from './modules/process/QC/SampleRunMultiQC'
include { QcPileup }                           from './modules/process/QC/QcPileup'
include { QcConpair }                          from './modules/process/QC/QcConpair'
include { SomaticRunMultiQC }                  from './modules/process/QC/SomaticRunMultiQC'
include { QcConpairAll }                       from './modules/process/QC/QcConpairAll'
include { CohortRunMultiQC }                   from './modules/process/QC/CohortRunMultiQC'

//Aggregate
include { GermlineAggregateMaf }               from './modules/process/Aggregate/GermlineAggregateMaf'
include { GermlineAggregateSv }                from './modules/process/Aggregate/GermlineAggregateSv'
include { QcBamAggregate }                     from './modules/process/Aggregate/QcBamAggregate'
include { QcConpairAggregate }                 from './modules/process/Aggregate/QcConpairAggregate'
include { SomaticAggregateFacets }             from './modules/process/Aggregate/SomaticAggregateFacets'
include { SomaticAggregateLOHHLA }             from './modules/process/Aggregate/SomaticAggregateLOHHLA'
include { SomaticAggregateMaf }                from './modules/process/Aggregate/SomaticAggregateMaf'
include { SomaticAggregateMetadata }           from './modules/process/Aggregate/SomaticAggregateMetadata'
include { SomaticAggregateNetMHC }             from './modules/process/Aggregate/SomaticAggregateNetMHC'
include { SomaticAggregateSv }                 from './modules/process/Aggregate/SomaticAggregateSv'


println ''

pairingQc = params.pairing

referenceMap = defineReferenceMap()
targetsMap   = loadTargetReferences()

workflow alignment_wf
{
  main:
    if (params.bamMapping)
    {
      println "Alignment workflow cannot accept bam files for input."
      exit 1
    }
    if(params.mapping)
    {
      TempoUtils.checkAssayType(params.assayType)
      if (params.watch == false) {
        keySet = targetsMap.keySet()
        mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)
        inputMapping = params.mapping ? TempoUtils.extractFastq(mappingFile, params.assayType, targetsMap.keySet()) : TempoUtils.extractBAM(mappingFile, params.assayType, targetsMap.keySet())
      }
      else if (params.watch == true) {
        mappingFile = params.mapping ? file(params.mapping, checkIfExists: false) : file(params.bamMapping, checkIfExists: false)
        inputMapping  = params.mapping ? watchMapping(mappingFile, params.assayType) : watchBamMapping(mappingFile, params.assayType)
        epochMap[params.mapping ? params.mapping : params.bamMapping ] = 0
      }
      else{}

/*      if(params.pairing){
        if (params.watch == false) {
          pairingFile = file(params.pairing, checkIfExists: true)
          inputPairing  = TempoUtils.extractPairing(pairingFile)
          TempoUtils.crossValidateTargets(inputMapping, inputPairing)

          samplesInMapping = inputMapping.flatMap{[it[0]]}.unique().toSortedList()
          samplesInPairing = inputPairing.flatten().unique().toSortedList()
          CrossValidateSamples(samplesInMapping, samplesInPairing)
          if(!CrossValidateSamples.out.isValid) {exit 1}
        }
        else if (params.watch == true) {
          pairingFile = file(params.pairing, checkIfExists: false)
          inputPairing  = watchPairing(pairingFile)
          epochMap[params.pairing] = 0 
        }
        else{}
      }
*/

      // Parse input FASTQ mapping
      if (params.watch != true) {
        inputMapping.groupTuple(by: [0])
                  .map { idSample, targets, files_pe1, files_pe2
                    -> tuple(groupKey(idSample, targets.size()), targets, files_pe1, files_pe2)
                  }
                  .transpose()
            .set { inputMapping }
      }

      inputMapping.map { idSample, target, file_pe1, file_pe2 ->
                      [idSample, target, file_pe1, file_pe2, idSample + '@' + file_pe1.getSimpleName(), file_pe1.getSimpleName()]
      }
                  .set { inputFastqs }

      if (params.splitLanes) {
        inputFastqs.set { fastqsNeedSplit }
        inputFastqs.set { fastqsNoNeedSplit }

        fastqsNeedSplit
              .filter { item -> !(item[2].getName() =~ /_L(\d){3}_/) }
              .multiMap { idSample, target, file_pe1, file_pe2, fileID, lane ->
                inputFastqR1: [idSample, target, file_pe1, file_pe1.toString()]
                inputFastqR2: [idSample, target, file_pe2, file_pe2.toString()]
              }
              .set { fastqsNeedSplit }

        fastqsNoNeedSplit
              .filter { item -> item[2].getName() =~ /_L(\d){3}_/ }
              .map { idSample, target, file_pe1, file_pe2, fileID, lane
                -> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
              }
              .set { fastqsNoNeedSplit }

        perLaneFastqsR1 = SplitLanesR1(fastqsNeedSplit.inputFastqR1).R1SplitData
        perLaneFastqsR2 = SplitLanesR2(fastqsNeedSplit.inputFastqR2).R2SplitData

        def fastqR1fileIDs = [:]
        perLaneFastqsR1 = perLaneFastqsR1.transpose()
              .map { item ->
                  def idSample = item[0]
                  def target = item[1]
                  def fastq = item[2]
                  def fileID = idSample + '@' + item[3].getSimpleName()
                  def lane = fastq.getSimpleName().split('_L00')[1].split('_')[0]
                  def laneCount = item[4].getSimpleName().toInteger()

                  // This only checks if same read groups appears in two or more fastq files which belongs to the same sample. Cross sample check will be performed after AlignReads since the read group info is not available for fastqs which does not need to be split.
                  if ( !params.watch ) {
                    if (!TempoUtils.checkDuplicates(fastqR1fileIDs, fileID + '@' + lane, fileID + "\t" + fastq, 'the following fastq files since they contain the same RGID')) { exit 1 }
                  }
                  [idSample, target, fastq, fileID, lane, laneCount]
              }

        def fastqR2fileIDs = [:]
        perLaneFastqsR2 = perLaneFastqsR2.transpose()
                  .map { item ->
                    def idSample = item[0]
                    def target = item[1]
                    def fastq = item[2]
                    def fileID = idSample + '@' + item[3].getSimpleName()
                    def lane = fastq.getSimpleName().split('_L00')[1].split('_')[0]
                    def laneCount = item[4].getSimpleName().toInteger()
                    if ( !params.watch ) {
                      if (!TempoUtils.checkDuplicates(fastqR2fileIDs, fileID + '@' + lane, fileID + "\t" + fastq, 'the follwoing fastq files since they contain the same RGID')) { exit 1 }
                    }
                    [idSample, target, fastq, fileID, lane, laneCount]
                  }

        fastqFiles  = perLaneFastqsR1
              .mix(perLaneFastqsR2)
              .groupTuple(by: [0, 1, 3, 4, 5], size: 2, sort: true)
              .map {  idSample, target, fastqPairs, fileID, lanes, laneCount ->
                tuple(idSample, target, fastqPairs, groupKey(fileID, laneCount), lanes)
              }
              .map { idSample, target, fastqPairs, fileID, lane ->
                  [idSample, target, fastqPairs[0], fastqPairs[1], fileID, lane]
              }
              .map { item ->
                def idSample = item[0]
                def target = item[1]
                def fastqPair1 = item[2]
                def fastqPair2 = item[3]
                if (item[2].toString().split('_R1').size() < item[3].toString().split('_R1').size()) {
                  fastqPair1 = item[3]
                  fastqPair2 = item[2]
                }
                def fileID = item[4]
                def lane = item[5]
                [idSample, target, fastqPair1, fastqPair1.size(), fastqPair2, fastqPair2.size(), fileID, lane]
              }
              .mix(fastqsNoNeedSplit)
      }
      else {
        fastqFiles = inputFastqs.map { idSample, target, file_pe1, file_pe2, fileID, lane
          -> tuple(idSample, target, file_pe1, file_pe1.size(),
                                              file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
        }
      }

      //Align reads to reference.
      AlignReads(fastqFiles, Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex]))

      AlignReads.out.fastPJson4MultiQC
        .groupTuple(by:[2])
        .map{idSample, jsonFile, fileID -> 
          def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
          [idSampleout, jsonFile]
        }.groupTuple(by: [0])
        .map{ idSample, jsonFile -> 
          [idSample, jsonFile.flatten()]
        }.set{ fastPJson } 

      // Check for FASTQ files which might have different path but contains the same reads, based only on the name of the first read.
      def allReadIds = [:]
      AlignReads.out.sortedBam.map { idSample, target, bam, fileID, lane, readIdFile -> def readId = '@' + readIdFile.getSimpleName().replaceAll('@', ':')
        // Use the first line of the fastq file (the name of the first read) as unique identifier to check across all the samples if there is any two fastq files contains the same read name, if so, we consider there are some human error of mixing up the same reads into different fastq files
        if ( !params.watch ) {
          if (!TempoUtils.checkDuplicates(allReadIds, readId, idSample + "\t" + bam, "the following samples, since they contain the same read: \n${ readId }")) {exit 1}
        }
        [idSample, target, bam, fileID, lane]
      }
      .groupTuple(by: [3])
      .map { item ->
        def idSample = item[0] instanceof Collection ? item[0].first() : item[0]
        def target   = item[1] instanceof Collection ? item[1].first() : item[1]
        def bams = item[2]
        [idSample, target, bams]
      }
      .groupTuple(by: [0])
      .map { item ->
        def idSample = item[0]
        def target =  item[1] instanceof Collection ? item[1].first() : item[1]
        def bams = item[2].flatten()
        [idSample, bams, target]
      }
      .set { groupedBam }

      MergeBamsAndMarkDuplicates(groupedBam)
      RunBQSR(MergeBamsAndMarkDuplicates.out.mdBams,
              Channel.value([
                referenceMap.genomeFile,
                referenceMap.genomeIndex,
                referenceMap.genomeDict,
                referenceMap.dbsnp,
                referenceMap.dbsnpIndex,
                referenceMap.knownIndels,
                referenceMap.knownIndelsIndex
              ]))


      File file = new File(outname)
      file.newWriter().withWriter { w ->
          w << "SAMPLE\tTARGET\tBAM\tBAI\n"
      }

      RunBQSR.out.bamResults.subscribe { Object obj ->
        file.withWriterAppend { out ->
            out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}"
        }
      }
    }
    else{
      if(params.pairing){
        println "ERROR: When --pairing [tsv], --mapping [tsv] must be provided."
        exit 1
      }
    }
  

  emit:
    //CrossValidateSamples_output = CrossValidateSamples.out
    //SplitLanesR1_output = SplitLanesR1.out
    //SplitLanesR2_output = SplitLanesR2.out
    //AlignReads_output = AlignReads.out
    //MergeBamsAndMarkDuplicates_output = MergeBamsAndMarkDuplicates.out
    RunBQSR_bamsBQSR = RunBQSR.out.bamsBQSR
    RunBQSR_bamResults = RunBQSR.out.bamResults
    RunBQSR_bamSize = RunBQSR.out.bamSize

}

workflow {
  if (params.mapping || params.bamMapping) {
    TempoUtils.checkAssayType(params.assayType)
    if (params.watch == false) {
      keySet = targetsMap.keySet()
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)
      inputMapping = params.mapping ? TempoUtils.extractFastq(mappingFile, params.assayType, targetsMap.keySet()) : TempoUtils.extractBAM(mappingFile, params.assayType, targetsMap.keySet())
    }
    else if (params.watch == true) {
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: false) : file(params.bamMapping, checkIfExists: false)
      inputMapping  = params.mapping ? watchMapping(mappingFile, params.assayType) : watchBamMapping(mappingFile, params.assayType)
      epochMap[params.mapping ? params.mapping : params.bamMapping ] = 0
    }
    else{}
    if(params.pairing){
      if(runQC){
        pairingQc = true
      }
      if (params.watch == false) {
        pairingFile = file(params.pairing, checkIfExists: true)
        inputPairing  = TempoUtils.extractPairing(pairingFile)
        TempoUtils.crossValidateTargets(inputMapping, inputPairing)

        samplesInMapping = inputMapping.flatMap{[it[0]]}.unique().toSortedList()
        samplesInPairing = inputPairing.flatten().unique().toSortedList()
        CrossValidateSamples(samplesInMapping, samplesInPairing)
        if(!CrossValidateSamples.out.isValid) {exit 1}
      }
      else if (params.watch == true) {
        pairingFile = file(params.pairing, checkIfExists: false)
        inputPairing  = watchPairing(pairingFile)
        epochMap[params.pairing] = 0 
      }
      else{}
      /*if (!runSomatic && !runGermline && !runQC){
        println "ERROR: --pairing [tsv] is not used because none of --somatic/--germline/--QC was enabled. If you only need to do BAM QC and/or BAM generation, remove --pairing [tsv]."
        exit 1
      }
      */
    }
    else{
      if (runSomatic || runGermline){
        println "ERROR: --pairing [tsv] needed when using --mapping/--bamMapping [tsv] with --somatic/--germline"
        exit 1
      }
    }
  }
  else{
    if(params.pairing){
      println "ERROR: When --pairing [tsv], --mapping/--bamMapping [tsv] must be provided."
      exit 1
    }
  }

  if (!runSomatic && runGermline) {
    println "WARNING: You can't run GERMLINE section without running SOMATIC section. Activating SOMATIC section automatically"
    runSomatic = true
  }
  if (runAggregate == false) {
    if (!params.mapping && !params.bamMapping) {
      println 'ERROR: (--mapping/-bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv] ) or (--aggregate [tsv]) need to be provided, otherwise nothing to be run.'
      exit 1
    }
  }
  else if (runAggregate == true) {
    if ((params.mapping || params.bamMapping) && params.pairing) {
      if (!(runSomatic || runGermline || runQC)) {
        println 'ERROR: Nothing to be aggregated. One or more of the option --somatic/--germline/--QC need to be enabled when using --aggregate'
      }
    }
    else if ((params.mapping || params.bamMapping) && !params.pairing) {
      if (!runQC) {
        println 'ERROR: Nothing to be aggregated. --QC need to be enabled when using --mapping/--bamMapping [tsv], --pairing false and --aggregate true.'
        exit 1
      }
    }
    else {
      println 'ERROR: (--mapping/--bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv]) or (--aggregate [tsv]) need to be provided when using --aggregate true'
      println '       If you want to run aggregate only, you need to use --aggregate [tsv]. See manual'
      exit 1
    }
  }
  else {
    if ((runSomatic || runGermline || runQC) && !params.mapping && !params.bamMapping) {
      println 'ERROR: Conflict input! When running --aggregate [tsv] with --mapping/--bamMapping/--pairing [tsv] disabled, --QC/--somatic/--germline all need to be disabled!'
      println "       If you want to run aggregate somatic/germline/qc, just include an additianl colum PATH in the [tsv] and no need to use --QC/--somatic/--germline flag, since it's auto detected. See manual"
      exit 1
    }
  }

  if (!(params.cosmic in ['v2', 'v3'])) {
    println "ERROR: Possible values of mutational signature reference --cosmic is 'v2', 'v3'"
    exit 1
  }


  // Skip these processes if starting from aligned BAM files
  if (params.alignWF)
  {
    alignment_wf()
  }

 /* if (params.mapping) {
    // Parse input FASTQ mapping
    if (params.watch != true) {
      inputMapping.groupTuple(by: [0])
                .map { idSample, targets, files_pe1, files_pe2
                  -> tuple(groupKey(idSample, targets.size()), targets, files_pe1, files_pe2)
                }
                .transpose()
          .set { inputMapping }
    }

    inputMapping.map { idSample, target, file_pe1, file_pe2 ->
                    [idSample, target, file_pe1, file_pe2, idSample + '@' + file_pe1.getSimpleName(), file_pe1.getSimpleName()]
    }
                .set { inputFastqs }

    if (params.splitLanes) {
      inputFastqs.set { fastqsNeedSplit }
      inputFastqs.set { fastqsNoNeedSplit }

      fastqsNeedSplit
            .filter { item -> !(item[2].getName() =~ /_L(\d){3}_/) }
            .multiMap { idSample, target, file_pe1, file_pe2, fileID, lane ->
              inputFastqR1: [idSample, target, file_pe1, file_pe1.toString()]
              inputFastqR2: [idSample, target, file_pe2, file_pe2.toString()]
            }
            .set { fastqsNeedSplit }

      fastqsNoNeedSplit
            .filter { item -> item[2].getName() =~ /_L(\d){3}_/ }
            .map { idSample, target, file_pe1, file_pe2, fileID, lane
              -> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
            }
            .set { fastqsNoNeedSplit }

      perLaneFastqsR1 = SplitLanesR1(fastqsNeedSplit.inputFastqR1).R1SplitData
      perLaneFastqsR2 = SplitLanesR2(fastqsNeedSplit.inputFastqR2).R2SplitData

      def fastqR1fileIDs = [:]
      perLaneFastqsR1 = perLaneFastqsR1.transpose()
            .map { item ->
                def idSample = item[0]
                def target = item[1]
                def fastq = item[2]
                def fileID = idSample + '@' + item[3].getSimpleName()
                def lane = fastq.getSimpleName().split('_L00')[1].split('_')[0]
                def laneCount = item[4].getSimpleName().toInteger()

                // This only checks if same read groups appears in two or more fastq files which belongs to the same sample. Cross sample check will be performed after AlignReads since the read group info is not available for fastqs which does not need to be split.
                if ( !params.watch ) {
                  if (!TempoUtils.checkDuplicates(fastqR1fileIDs, fileID + '@' + lane, fileID + "\t" + fastq, 'the following fastq files since they contain the same RGID')) { exit 1 }
                }
                [idSample, target, fastq, fileID, lane, laneCount]
            }

      def fastqR2fileIDs = [:]
      perLaneFastqsR2 = perLaneFastqsR2.transpose()
                .map { item ->
                  def idSample = item[0]
                  def target = item[1]
                  def fastq = item[2]
                  def fileID = idSample + '@' + item[3].getSimpleName()
                  def lane = fastq.getSimpleName().split('_L00')[1].split('_')[0]
                  def laneCount = item[4].getSimpleName().toInteger()
                  if ( !params.watch ) {
                    if (!TempoUtils.checkDuplicates(fastqR2fileIDs, fileID + '@' + lane, fileID + "\t" + fastq, 'the follwoing fastq files since they contain the same RGID')) { exit 1 }
                  }
                  [idSample, target, fastq, fileID, lane, laneCount]
                }

      fastqFiles  = perLaneFastqsR1
            .mix(perLaneFastqsR2)
            .groupTuple(by: [0, 1, 3, 4, 5], size: 2, sort: true)
            .map {  idSample, target, fastqPairs, fileID, lanes, laneCount ->
              tuple(idSample, target, fastqPairs, groupKey(fileID, laneCount), lanes)
            }
            .map { idSample, target, fastqPairs, fileID, lane ->
                [idSample, target, fastqPairs[0], fastqPairs[1], fileID, lane]
            }
            .map { item ->
              def idSample = item[0]
              def target = item[1]
              def fastqPair1 = item[2]
              def fastqPair2 = item[3]
              if (item[2].toString().split('_R1').size() < item[3].toString().split('_R1').size()) {
                fastqPair1 = item[3]
                fastqPair2 = item[2]
              }
              def fileID = item[4]
              def lane = item[5]
              [idSample, target, fastqPair1, fastqPair1.size(), fastqPair2, fastqPair2.size(), fileID, lane]
            }
            .mix(fastqsNoNeedSplit)
    }
    else {
      fastqFiles = inputFastqs.map { idSample, target, file_pe1, file_pe2, fileID, lane
        -> tuple(idSample, target, file_pe1, file_pe1.size(),
                                             file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
      }
    }

    //Align reads to reference.
    AlignReads(fastqFiles, Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex]))

    AlignReads.out.fastPJson4MultiQC
      .groupTuple(by:[2])
      .map{idSample, jsonFile, fileID -> 
        def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
        [idSampleout, jsonFile]
      }.groupTuple(by: [0])
      .map{ idSample, jsonFile -> 
        [idSample, jsonFile.flatten()]
      }.set{ fastPJson } 

    // Check for FASTQ files which might have different path but contains the same reads, based only on the name of the first read.
    def allReadIds = [:]
    AlignReads.out.sortedBam.map { idSample, target, bam, fileID, lane, readIdFile -> def readId = '@' + readIdFile.getSimpleName().replaceAll('@', ':')
      // Use the first line of the fastq file (the name of the first read) as unique identifier to check across all the samples if there is any two fastq files contains the same read name, if so, we consider there are some human error of mixing up the same reads into different fastq files
      if ( !params.watch ) {
        if (!TempoUtils.checkDuplicates(allReadIds, readId, idSample + "\t" + bam, "the following samples, since they contain the same read: \n${ readId }")) {exit 1}
      }
      [idSample, target, bam, fileID, lane]
    }
    .groupTuple(by: [3])
    .map { item ->
      def idSample = item[0] instanceof Collection ? item[0].first() : item[0]
      def target   = item[1] instanceof Collection ? item[1].first() : item[1]
      def bams = item[2]
      [idSample, target, bams]
    }
    .groupTuple(by: [0])
    .map { item ->
      def idSample = item[0]
      def target =  item[1] instanceof Collection ? item[1].first() : item[1]
      def bams = item[2].flatten()
      [idSample, bams, target]
    }
    .set { groupedBam }

    MergeBamsAndMarkDuplicates(groupedBam)
    RunBQSR(MergeBamsAndMarkDuplicates.out.mdBams,
            Channel.value([
              referenceMap.genomeFile,
              referenceMap.genomeIndex,
              referenceMap.genomeDict,
              referenceMap.dbsnp,
              referenceMap.dbsnpIndex,
              referenceMap.knownIndels,
              referenceMap.knownIndelsIndex
            ]))


    File file = new File(outname)
    file.newWriter().withWriter { w ->
        w << "SAMPLE\tTARGET\tBAM\tBAI\n"
    }

    RunBQSR.out.bamResults.subscribe { Object obj ->
      file.withWriterAppend { out ->
          out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}"
      }
    }
  }
  */
   //end if params.mapping


  /*
  ================================================================================
  =                                PAIRING TUMOR and NORMAL                      =
  ================================================================================
  */
  // If starting with BAM files, parse BAM pairing input
 
  if (params.bamMapping) {
    inputChannel = inputMapping
    if (runQC){
      inputMapping.map{idSample, target, bam, bai ->
        [ idSample,target, bam.getParent() ]
      }.set{ locateFastP4MultiQC }

      locateFastP4MultiQC.map{ idSample,target, bamFolder -> 
        [idSample, file(bamFolder + "/fastp/*json")]
      }.set{ fastPJson }
    } 
  }
  else
  {
    //inputChannel = RunBQSR.out.bamsBQSR
    inputChannel = alignment_wf.out.RunBQSR_bamsBQSR
  }

  if (params.pairing) {
    // Parse input FASTQ mapping
    inputChannel.combine(inputPairing)
      .filter { item ->
        def idSample = item[0]
        def target = item[1]
        def sampleBam = item[2]
        def sampleBai = item[3]
        def idTumor = item[4]
        def idNormal = item[5]
        idSample == idTumor
      }.map { item ->
        def idTumor = item[4]
        def idNormal = item[5]
        def tumorBam = item[2]
        def tumorBai = item[3]
        def target = item[1]
        return [ idTumor, idNormal, target, tumorBam, tumorBai ]
      }
      .unique()
      .set{ bamsTumor }

    inputChannel.combine(inputPairing)
      .filter { item ->
        def idSample = item[0]
        def target = item[1]
        def sampleBam = item[2]
        def sampleBai = item[3]
        def idTumor = item[4]
        def idNormal = item[5]
        idSample == idNormal
      }.map { item ->
        def idTumor = item[4]
        def idNormal = item[5]
        def normalBam = item[2]
        def normalBai = item[3]
        def target = item[1]
        return [ idTumor, idNormal, target, normalBam, normalBai ]
      }
      .unique()
      .set{ bamsNormal }

    bamsNormal.map { item ->
      def idNormal = item[1]
      def target = item[2]
      def normalBam = item[3]
      def normalBai = item[4]
      return [ idNormal, target, normalBam, normalBai ] }
      .unique()
      .set{ bams }

    bamsTumor.combine(bamsNormal, by: [0,1,2])
      .map { item -> // re-order the elements
        def idTumor = item[0]
        def idNormal = item[1]
        def target = item[2]
        def bamTumor = item[3]
        def baiTumor = item[4]
        def bamNormal = item[5]
        def baiNormal = item[6]

        return [ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal ]
      }
      .set{ bamFiles }

  } // End of "if (pairingQc) {}"


  /*
  ================================================================================
  =                                SOMATIC PIPELINE                              =
  ================================================================================
  */
  if (runSomatic || runGermline || runQC) {
    // parse --tools parameter for downstream 'when' conditionals, e.g. when: 'delly ' in tools
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

    // If using running either conpair or conpairAll, run pileup as well to generate pileups
    if ("conpair" in tools) {
      tools.add("pileup")
    }
    if ("conpairall" in tools) {
      runConpairAll = true
      tools.add("pileup")
    }
  }

  if (runSomatic || runGermline) {
    targets4Intervals = Channel.from(targetsMap.keySet())
      .map{ targetId ->
        [ targetId, targetsMap[targetId].targetsBedGz, targetsMap[targetId].targetsBedGzTbi ]
      }

    targets4Intervals = Channel.from(targetsMap.keySet())
      .map{ targetId ->
        [ targetId, targetsMap."${targetId}".targetsBedGz, targetsMap."${targetId}".targetsBedGzTbi ]
      }

    CreateScatteredIntervals(Channel.value([referenceMap.genomeFile, 
                                            referenceMap.genomeIndex, 
                                            referenceMap.genomeDict]), 
                                            targets4Intervals,
                                            runSomatic, runGermline)

    //Associating interval_list files with BAM files, putting them into one channel
    bamFiles.combine(CreateScatteredIntervals.out.mergedIList, by: 2)
      .map{
        item ->
          def idTumor = item[1]
          def idNormal = item[2]
          def target = item[0]
          def tumorBam = item[3]
          def normalBam = item[4]
          def tumorBai = item[5]
          def normalBai = item[6]
          def intervalBed = item[7]
          def key = idTumor+"__"+idNormal+"@"+target // adding one unique key

          return [ key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed ]
      }.map{ 
        key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed -> 
        tuple ( 
            groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
            idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed
        )
    }
    .transpose()
    .set{ mergedChannelSomatic }

    bams.combine(CreateScatteredIntervals.out.mergedIList, by: 1)
      .map{
        item ->
          def idNormal = item[1]
          def target = item[0]
          def normalBam = item[2]
          def normalBai = item[3]
          def intervalBed = item[4]
          def key = idNormal+"@"+target // adding one unique key

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
  } //end if (runSomatic || runGermline)


  if (runSomatic){
    RunMutect2(mergedChannelSomatic, 
              Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
              tools,
              runSomatic)

    //Formatting the channel to be keyed by idTumor, idNormal, and target
    // group by groupKey(key, intervalBed.size())
    RunMutect2.out.forMutect2Combine.groupTuple().set{ forMutect2Combine }
    SomaticCombineMutect2Vcf(forMutect2Combine, 
                             Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                             tools,
                             runSomatic)

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }
    SomaticDellyCall(svTypes,
                     bamFiles,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions]),
                     tools,
                     runSomatic)


    SomaticRunManta(bamFiles, 
                  Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                  Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex]),
                  tools,
                  runSomatic)


    // Put manta output and delly output into the same channel so they can be processed together in the group key
    // that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
    SomaticDellyCall.out.dellyFilter4Combine.groupTuple(by: [0,1,2], size: 5).combine(SomaticRunManta.out.manta4Combine, by: [0,1,2]).set{ dellyMantaCombineChannel }

    // --- Process Delly and Manta VCFs 
    // Merge VCFs, Delly and Manta
    SomaticMergeDellyAndManta(dellyMantaCombineChannel,
                              tools,
                              runSomatic)

    bamFiles.combine(SomaticRunManta.out.mantaToStrelka, by: [0, 1, 2])
        .map{ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi ->
              [idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi]
        }.set{ input4Strelka }

    SomaticRunStrelka2(input4Strelka,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                        tools,
                        runSomatic)

    SomaticCombineMutect2Vcf.out.mutect2CombinedVcfOutput.combine(bamFiles, by: [0,1,2]).combine(SomaticRunStrelka2.out.strelka4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }

    SomaticCombineChannel(mutectStrelkaChannel,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.exomePoN, referenceMap.wgsPoN,referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]),
                          tools,
                          runSomatic)

    SomaticAnnotateMaf(SomaticCombineChannel.out.mutationMergedVcf,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
                                        referenceMap.vepCache, referenceMap.isoforms
                                      ]),
                        tools,
                        runSomatic
                      )

    RunMutationSignatures(SomaticAnnotateMaf.out.mafFile, 
                          tools,
                          runSomatic)


    DoFacets(bamFiles,
              Channel.value([referenceMap.facetsVcf]),
              tools,
              runSomatic)

    DoFacetsPreviewQC(DoFacets.out.Facets4FacetsPreview,
                      tools,
                      runSomatic)

    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4Aggregate } // idTumor, idNormal, summaryFiles, qcFiles
    DoFacets.out.FacetsRunSummary.combine(DoFacetsPreviewQC.out.FacetsPreviewOut, by:[0,1]).set{ FacetsQC4SomaticMultiQC } // idTumor, idNormal, summaryFiles, qcFiles
    FacetsQC4Aggregate.map{ idTumor, idNormal, summaryFiles, qcFiles ->
      ["placeholder",idTumor, idNormal, summaryFiles, qcFiles]
    }.set{FacetsQC4Aggregate}

    RunPolysolver(bams,
                  tools,
                  runSomatic)

    bamFiles.combine(DoFacets.out.facetsPurity, by: [0,1,2])
      .combine(RunPolysolver.out.hlaOutput, by: [1,2])
      .set{ mergedChannelLOHHLA }

    RunLOHHLA(mergedChannelLOHHLA, 
              Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat]),
              tools,
              runSomatic)

    RunPolysolver.out.hlaOutput.combine(SomaticAnnotateMaf.out.mafFile, by: [1,2]).set{ input4Neoantigen }

    RunNeoantigen(input4Neoantigen,
                  Channel.value([referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS]),
                  tools,
                  runSomatic)

    DoFacets.out.facetsForMafAnno.combine(RunNeoantigen.out.mafFileForMafAnno, by: [0,1,2]).set{ facetsMafFileSomatic }

    SomaticFacetsAnnotation(facetsMafFileSomatic,
                            tools,
                            runSomatic)

    RunMsiSensor(bamFiles,
                 Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict, referenceMap.msiSensorList]),
                 tools,
                 runSomatic)

    DoFacets.out.facetsPurity.combine(SomaticFacetsAnnotation.out.maf4MetaDataParser, by: [0,1,2])
			   .combine(DoFacets.out.FacetsQC4MetaDataParser, by: [0,1,2])
			   .combine(RunMsiSensor.out.msi4MetaDataParser, by: [0,1,2])
			   .combine(RunMutationSignatures.out.mutSig4MetaDataParser, by: [0,1,2])
			   .combine(RunPolysolver.out.hlaOutput, by: [1,2])
			   .unique()
         .map{ idNormal, target, idTumor, purityOut, mafFile, qcOutput, msifile, mutSig, placeHolder, polysolverFile ->
          [idNormal, target, idTumor, purityOut, mafFile, qcOutput, msifile, mutSig, placeHolder, polysolverFile, targetsMap."$target".codingBed]
         }.set{ mergedChannelMetaDataParser }

    MetaDataParser(mergedChannelMetaDataParser,
                   runSomatic)

  }
  else { 
    if (params.pairing) {
      inputPairing.map{ idTumor, idNormal -> 
        ["placeHolder",idTumor, idNormal,"",""]
      }.set{ FacetsQC4Aggregate }
    }
  } //End of 'if runSomatic'


/*
================================================================================
=                                GERMLINE PIPELINE                              =
================================================================================
*/
  if (runGermline){
    GermlineRunHaplotypecaller(mergedChannelGermline,
                              Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                              tools,
                              runGermline)

    GermlineRunHaplotypecaller.out.haplotypecaller4Combine.groupTuple().set{ haplotypecaller4Combine }

    GermlineCombineHaplotypecallerVcf(haplotypecaller4Combine,
                                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                                      tools,
                                      runGermline)

    bams.map{ idNormal, target, bamNormal, baiNormal -> 
              [idNormal, target, bamNormal, baiNormal, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi]
            }.set{bamsForStrelkaGermline}

    GermlineRunStrelka2(bamsForStrelkaGermline, 
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                        tools,
                        runGermline)

    // Join HaploTypeCaller and Strelka outputs, bcftools.
    GermlineCombineHaplotypecallerVcf.out.haplotypecallerCombinedVcfOutput.combine(GermlineRunStrelka2.out.strelkaOutputGermline, by: [0,1,2])
            .combine(bamsTumor, by: [1,2])
            .set{ mergedChannelVcfCombine }

    GermlineCombineChannel(mergedChannelVcfCombine,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex,]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex, referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]),
                          tools,
                          runGermline)

    GermlineAnnotateMaf(GermlineCombineChannel.out.mutationMergedGermline,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
                                       referenceMap.vepCache, referenceMap.isoforms]),
                                       tools,
                                       runGermline)
    
    DoFacets.out.facetsForMafAnno.combine(GermlineAnnotateMaf.out.mafFileGermline, by: [0,1,2])
            .set{ facetsMafFileGermline }

    GermlineFacetsAnnotation(facetsMafFileGermline,
                             tools,
                             runGermline)

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypesGermline }

    GermlineDellyCall(svTypesGermline,
                      bams,
                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions]),
                      tools,
                      runGermline)

    GermlineRunManta(bams,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                     Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex]),
                     tools,
                     runGermline)

    GermlineDellyCall.out.dellyFilter4CombineGermline.groupTuple(by: [0,1], size: 5)
            .combine(GermlineRunManta.out.mantaOutputGermline, by: [0,1])
            .set{ dellyMantaChannelGermline }

    GermlineMergeDellyAndManta(dellyMantaChannelGermline,
                               Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                               tools,
                               runGermline)

  } //End of 'if runGermline'


/*
================================================================================
=                              Quality Control                                 =
================================================================================
*/

  if (runQC) {

    inputChannel.map{ idSample, target, bam, bai ->
        [idSample, target, bam, bai, targetsMap."$target".targetsInterval,  targetsMap."$target".baitsInterval]
    }.set{ bamsBQSR4HsMetrics }

    QcCollectHsMetrics(bamsBQSR4HsMetrics,
                       Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                       params.assayType,
                       runQC)

    if (runQC && params.assayType != "exome"){
      QcCollectHsMetrics.out.collectHsMetricsOutput
        .map{ idSample, target, bam, bai, targetList, baitList -> [idSample, ""]}
        .set{ collectHsMetricsOutput }
    }

    inputChannel
      .map{ idSample, target, bam, bai -> [ idSample, target, bam, bai, file(targetsMap."$target".targetsBed) ]}
      .set{ bamsBQSR4Qualimap }

    QcQualimap(bamsBQSR4Qualimap,
               runQC)

    Channel.from(true, false).set{ ignore_read_groups }
    inputChannel
      .map{ idSample, target, bam, bai -> 
        [ idSample, target, bam, bai, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi ]
      }.set{ bamsBQSR4Alfred }

    QcAlfred(ignore_read_groups, 
             bamsBQSR4Alfred,
             Channel.value([referenceMap.genomeFile]),
             runQC)

    QcAlfred.out.alfredOutput
      .groupTuple(size:2, by:0)
      .join(fastPJson, by:0)
      .join(QcQualimap.out.qualimap4Process, by:0)
      .join(QcCollectHsMetrics.out.collectHsMetricsOutput, by:0)
      .set{ sampleMetrics4MultiQC }

    SampleRunMultiQC(sampleMetrics4MultiQC, 
                     Channel.value([multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo]))


    if (pairingQc) {
      QcPileup(inputChannel,
               Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
               tools,
               runQC)

      QcPileup.out.pileupOutput.combine(inputPairing)
              .filter { item ->
                def idSample = item[0]
                def samplePileup = item[1]
                def idTumor = item[2]
                def idNormal = item[3]
                idSample == idTumor
              }.map { item ->
                def idTumor = item[2]
                def idNormal = item[3]
                def tumorPileup = item[1]
                return [ idTumor, idNormal, tumorPileup ]
              }
              .unique()
              .set{ pileupT }

      QcPileup.out.pileupOutput.combine(inputPairing)
              .filter { item ->
                def idSample = item[0]
                def samplePileup = item[1]
                def idTumor = item[2]
                def idNormal = item[3]
                idSample == idNormal
              }.map { item ->
                def idTumor = item[2]
                def idNormal = item[3]
                def normalPileup = item[1]
                return [ idTumor, idNormal, normalPileup ]
              }
              .unique()
              .set{ pileupN }

      pileupT.combine(pileupN, by: [0, 1]).unique().set{ pileupConpair }

      QcConpair(pileupConpair,
                Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                tools,
                runQC)

      if (runSomatic) {
        QcConpair.out.conpairOutput
          .map{ placeHolder, idTumor, idNormal, conpairFiles -> [idTumor, idNormal, conpairFiles]}
          .join(FacetsQC4SomaticMultiQC, by:[0,1])
          .set{ somaticMultiQCinput }
      } else {
        QcConpair.out.conpairOutput
          .map{ placeHolder, idTumor, idNormal, conpairFiles -> [idTumor, idNormal, conpairFiles, "", ""]}
          .set{ somaticMultiQCinput }
      }

      SomaticRunMultiQC(somaticMultiQCinput,
                        Channel.value([multiqcWesConfig,multiqcWgsConfig,multiqcTempoLogo]))

      if(runConpairAll){
        pileupT.combine(pileupN).unique().set{ pileupConpairAll }

        QcConpairAll(pileupConpairAll,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                     runConpairAll,
                     runQC)
      }

      // -- Run based on QcConpairAll channels or the single QcConpair channels
      conpairConcord4Aggregate = (!runConpairAll ? QcConpair.out.conpairConcord : QcConpairAll.out.conpairAllConcord)
      conpairContami4Aggregate = (!runConpairAll ? QcConpair.out.conpairContami : QcConpairAll.out.conpairAllContami)

    } // End of "if (pairingQc)". Doing QcPileup or QcConpair/QcConpairAll only when --pairing [tsv] is given

  } // End of "if (runQc)" 




/*
================================================================================
=                              Cohort Aggregation                              =
================================================================================
*/
  if ( !params.mapping && !params.bamMapping ){
    runSomatic = true
    runGermline = true
    runQC = true
    pairingQc = true
    if(!params.watch){
      TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
        .set{ inputAggregate }
    }
    else{
      watchAggregateWithPath(file(runAggregate, checkIfExists: true))
        .set{ inputAggregate }
      epochMap[runAggregate] = 0
    }

    inputAggregate.multiMap{ cohort, idTumor, idNormal, path ->
		  finalMaf4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.final.maf" )]
      NetMhcStats4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.all_neoantigen_predictions.txt")]
      FacetsPurity4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*_purity.seg")]
      FacetsHisens4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*_hisens.seg")]
      FacetsOutLog4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*_OUT.txt")]
      FacetsArmLev4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*.arm_level.txt")]
      FacetsGeneLev4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*/*.gene_level.txt")]
      FacetsQC4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*_OUT.txt"), file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*/*.facets_qc.txt")]
      dellyMantaCombined4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.delly.manta.vcf.gz")]
      dellyMantaCombinedTbi4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.delly.manta.vcf.gz.tbi")]
      predictHLA4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.DNA.HLAlossPrediction_CI.txt")]
      intCPN4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*DNA.IntegerCPN_CI.txt")]
      MetaData4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.sample_data.txt")]
      mafFile4AggregateGermline: [idTumor, idNormal, cohort, "placeHolder", file(path + "/germline/" + idNormal + "/*/" + idTumor + "__" + idNormal + ".germline.final.maf")]
      dellyMantaCombined4AggregateGermline: [idNormal, cohort, idTumor, "placeHolder", "noTumor", file(path + "/germline/" + idNormal + "/*/*.delly.manta.vcf.gz")]
      dellyMantaCombinedTbi4AggregateGermline: [idNormal, cohort, idTumor, "placeHolder", "noTumor", file(path + "/germline/" + idNormal + "/*/*.delly.manta.vcf.gz.tbi")]
      conpairConcord4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/" + idTumor + "__" + idNormal + ".concordance.txt")]
      conpairContami4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/" + idTumor + "__" + idNormal + ".contamination.txt")]
		  alfredIgnoreYTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.alfred.tsv.gz/")]
		  alfredIgnoreYNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.alfred.tsv.gz/")]
		  alfredIgnoreNTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.alfred.per_readgroup.tsv.gz/")]
		  alfredIgnoreNNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.alfred.per_readgroup.tsv.gz/")]
      qualimapTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/qualimap/${idTumor}_qualimap_rawdata.tar.gz")]
      qualimapNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/qualimap/${idNormal}_qualimap_rawdata.tar.gz")]
		  hsMetricsTumor: params.assayType == "exome" ? [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.hs_metrics.txt")] : [cohort, idTumor, idNormal, ""]
		  hsMetricsNormal: params.assayType == "exome" ? [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.hs_metrics.txt")] : [cohort, idTumor, idNormal, ""]
      fastpTumor: [cohort, idTumor, idNormal, file(path + "/bams/" + idTumor + "/*/*.fastp.json")]
      fastpNormal: [cohort, idTumor, idNormal, file(path + "/bams/" + idNormal + "/*/*.fastp.json")]
    }
		.set { aggregateList }

    inputSomaticAggregateMaf = aggregateList.finalMaf4Aggregate.transpose().groupTuple(by:[2])
    inputSomaticAggregateNetMHC = aggregateList.NetMhcStats4Aggregate.transpose().groupTuple(by:[2])
    inputPurity4Aggregate = aggregateList.FacetsPurity4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputHisens4Aggregate = aggregateList.FacetsHisens4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputOutLog4Aggregate = aggregateList.FacetsOutLog4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputArmLev4Aggregate = aggregateList.FacetsArmLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputGeneLev4Aggregate = aggregateList.FacetsGeneLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputFacetsQC4CohortMultiQC = aggregateList.FacetsQC4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    inputSomaticAggregateSv = aggregateList.dellyMantaCombined4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputPredictHLA4Aggregate = aggregateList.predictHLA4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputIntCPN4Aggregate = aggregateList.intCPN4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateMetadata = aggregateList.MetaData4Aggregate.transpose().groupTuple(by:[2])
    inputGermlineAggregateMaf = aggregateList.mafFile4AggregateGermline.transpose().groupTuple(by:[2])
    inputGermlineAggregateSv = aggregateList.dellyMantaCombined4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputGermlineAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputAlfredIgnoreY = aggregateList.alfredIgnoreYTumor.unique().combine(aggregateList.alfredIgnoreYNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputAlfredIgnoreN = aggregateList.alfredIgnoreNTumor.unique().combine(aggregateList.alfredIgnoreNNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputQualimap4CohortMultiQC = aggregateList.qualimapTumor.unique().combine(aggregateList.qualimapNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputHsMetrics = aggregateList.hsMetricsTumor.unique().combine(aggregateList.hsMetricsNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    aggregateList.conpairConcord4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairConcord}
    aggregateList.conpairContami4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairContami}
    aggregateList.fastpTumor.unique().combine(aggregateList.fastpNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}.set{inputFastP4MultiQC}
  
  }
  else if(!(runAggregate == false)) {
    if (!(runAggregate == true)){
      if(!params.watch){
        TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
            .groupTuple()
            .map{ cohort, idTumor, idNormal, pathNoUse
                  -> tuple( groupKey(cohort, idTumor instanceof Collection ? idTumor.size() : 1), idTumor, idNormal)
            }
            .transpose()
            .set{ inputAggregate }
      }
      else{
        watchAggregate(file(runAggregate, checkIfExists: false))
          .set{ inputAggregate }
        epochMap[runAggregate] = 0
      }
    }
    else if(runAggregate == true){
      inputPairing.set{ cohortTable }
      cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
          .set{ inputAggregate }
    }
    else{}
              
    if (runSomatic){
      inputSomaticAggregateMaf = inputAggregate.combine(SomaticFacetsAnnotation.out.finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputSomaticAggregateNetMHC = inputAggregate.combine(RunNeoantigen.out.NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputPurity4Aggregate = inputAggregate.combine(DoFacets.out.FacetsPurity4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputHisens4Aggregate = inputAggregate.combine(DoFacets.out.FacetsHisens4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputOutLog4Aggregate = inputAggregate.combine(DoFacets.out.FacetsOutLog4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputArmLev4Aggregate = inputAggregate.combine(DoFacets.out.FacetsArmLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputGeneLev4Aggregate = inputAggregate.combine(DoFacets.out.FacetsGeneLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputSomaticAggregateSv = inputAggregate.combine(SomaticMergeDellyAndManta.out.dellyMantaCombined4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputSomaticAggregateSvTbi = inputAggregate.combine(SomaticMergeDellyAndManta.out.dellyMantaCombinedTbi4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputPredictHLA4Aggregate = inputAggregate.combine(RunLOHHLA.out.predictHLA4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputIntCPN4Aggregate = inputAggregate.combine(RunLOHHLA.out.intCPN4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputSomaticAggregateMetadata = inputAggregate.combine(MetaDataParser.out.MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])

      if (runGermline){
        inputGermlineAggregateMaf = inputAggregate.combine(GermlineFacetsAnnotation.out.mafFile4AggregateGermline, by:[1,2]).groupTuple(by:[2])
        inputGermlineAggregateSv = inputAggregate.combine(GermlineMergeDellyAndManta.out.dellyMantaCombined4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
        inputGermlineAggregateSvTbi = inputAggregate.combine(GermlineMergeDellyAndManta.out.dellyMantaCombinedTbi4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
      }
    }


    if (runQC){
      QcAlfred.out.bamsQcStats4Aggregate.branch{ item ->
            def idSample = item[0]
            def alfred = item[1]
      
            ignoreY: alfred =~ /.+\.alfred\.tsv\.gz/
            ignoreN: alfred =~ /.+\.alfred\.per_readgroup\.tsv\.gz/
        }
        .set{ bamsQcStats4Aggregate }

      inputPairing.combine(bamsQcStats4Aggregate.ignoreY)
            .branch { item ->
              def idTumor = item[0]
              def idNormal = item[1]
              def idSample = item[2]
              def alfred = item[3]
          
              tumor: idSample == idTumor
              normal: idSample == idNormal
            }
            .set{ alfredIgnoreY }

      alfredIgnoreY.tumor.combine(alfredIgnoreY.normal, by:[0,1])
            .combine(inputAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
            .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
            .groupTuple(by:[0])
            .map{ item ->
              def cohort = item[0]
              def idTumors = item[1].unique()
              def idNormals = item[2].unique()
              def fileTumor = item[3].unique()
              def fileNormal = item[4].unique()
        
              [cohort, fileTumor, fileNormal]
            }
            .unique()
            .set{ alfredIgnoreY }
      
      inputPairing.combine(bamsQcStats4Aggregate.ignoreN)
            .branch { item ->
              def idTumor = item[0]
              def idNormal = item[1]
              def idSample = item[2]
              def alfred = item[3]
          
              tumor: idSample == idTumor
              normal: idSample == idNormal
            }
            .set{ alfredIgnoreN }

      alfredIgnoreN.tumor.combine(alfredIgnoreN.normal, by:[0,1])
            .combine(inputAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
            .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
            .groupTuple(by:[0])
            .map{ item ->
              def cohort = item[0]
              def idTumors = item[1].unique()
              def idNormals = item[2].unique()
              def fileTumor = item[3].unique()
              def fileNormal = item[4].unique()
        
              [cohort, fileTumor, fileNormal]
            }
            .unique()
            .set{ alfredIgnoreN }
      
      inputPairing.combine(QcCollectHsMetrics.out.collectHsMetricsOutput)
            .branch { item ->
              def idTumor = item[0]
              def idNormal = item[1]
              def idSample = item[2]
              def hsMetrics = item[3]
          
              tumor: idSample == idTumor
              normal: idSample == idNormal
            }
            .set{ hsMetrics }

      hsMetrics.tumor.combine(hsMetrics.normal, by:[0,1])
              .combine(inputAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
              .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
              .groupTuple(by:[0])
              .map{ item ->
                def cohort = item[0]
                def idTumors = item[1].unique()
                def idNormals = item[2].unique()
                def fileTumor = item[3].unique()
                def fileNormal = item[4].unique()
                
                [cohort, fileTumor, fileNormal]
              }
              .unique()
              .set{ hsMetrics }

      inputHsMetrics = hsMetrics
      inputPairing.combine(fastPJson)
            .branch { item ->
              def idTumor = item[0]
              def idNormal = item[1]
              def idSample = item[2]
              def jsonFiles = item[3]
          
              tumor: idSample == idTumor
              normal: idSample == idNormal
            }
            .set{ fastPMetrics }

      fastPMetrics.tumor.combine(fastPMetrics.normal, by:[0,1])
              .combine(inputAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
              .map{ item -> [item[6], item[0], item[1], item[3], item[5]]}
              .groupTuple(by:[0])
              .map{ item ->
                def cohort = item[0]
                def idTumors = item[1].unique()
                def idNormals = item[2].unique()
                def fileTumor = item[3].flatten().unique()
                def fileNormal = item[4].flatten().unique()
                
                [cohort, fileTumor, fileNormal]
              }
              .unique()
              .set{ fastPMetrics }

      inputPairing.combine(QcQualimap.out.qualimap4Process)
        .branch { idTumor, idNormal, idSample, qualimapDir ->
          tumor:  idSample == idTumor
          normal: idSample == idNormal
        }
        .set{ qualimap4AggregateTN }

      qualimap4AggregateTN.tumor.combine(qualimap4AggregateTN.normal, by:[0,1])
        .combine(inputAggregate.map{ item -> [item[1], item[2], item[0]]}, by:[0,1])
        .map{ item -> [item[6], item[3], item[5]]}
        .groupTuple(by:[0])
        .map{ cohort, fileTumor, fileNormal ->
            [cohort, fileTumor.unique(), fileNormal.unique()]
        }
        .unique()
        .set{ inputQualimap4CohortMultiQC }

      inputAlfredIgnoreY = alfredIgnoreY
      inputAlfredIgnoreN = alfredIgnoreN
      inputFastP4MultiQC = fastPMetrics

      if (pairingQc){
        inputAggregate.combine(conpairConcord4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairConcord4Aggregate }
        inputAggregate.combine(conpairContami4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairContami4Aggregate }
        inputAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}.set{ inputFacetsQC4CohortMultiQC }
      }
    } //end 'if (runQC)'
  }//end 'else if(!(runAggregate == false))'
  else{}


  if (runAggregate && runSomatic) {
    SomaticAggregateMaf(inputSomaticAggregateMaf,
                        runSomatic)

    SomaticAggregateNetMHC(inputSomaticAggregateNetMHC,
                           runSomatic)

    inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
		     .join(inputOutLog4Aggregate, by:[0])
		     .join(inputArmLev4Aggregate, by:[0])
		     .join(inputGeneLev4Aggregate, by:[0])
		     .set{ inputSomaticAggregateFacets }

    SomaticAggregateFacets(inputSomaticAggregateFacets, 
                           runSomatic)

    inputSomaticAggregateSv.join(inputSomaticAggregateSvTbi)
		      .set{ inputSomaticAggregateSv }

    SomaticAggregateSv(inputSomaticAggregateSv,
                       runSomatic)

    inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
          .set{ inputSomaticAggregateLOHHLA }

    SomaticAggregateLOHHLA(inputSomaticAggregateLOHHLA,
                           runSomatic)

    SomaticAggregateMetadata(inputSomaticAggregateMetadata,
                             runSomatic)

  } //end 'if (runAggregate && runSomatic)'


  if (runAggregate && runGermline) {
    GermlineAggregateMaf(inputGermlineAggregateMaf,
                         runGermline)

    // --- Aggregate per-sample germline data, SVs
    inputGermlineAggregateSv.join(inputGermlineAggregateSvTbi)
          .set{ inputGermlineAggregateSv }

    GermlineAggregateSv(inputGermlineAggregateSv,
                        runGermline)
                        
  } 


  if (runAggregate && runQC) {
    inputAlfredIgnoreY.join(inputAlfredIgnoreN)
          .join(inputHsMetrics)
          .set{ inputQcBamAggregate }

    QcBamAggregate(inputQcBamAggregate,
                   runQC)

    if (pairingQc){
      inputConpairConcord4Aggregate.join(inputConpairContami4Aggregate)
			     .set{ inputQcConpairAggregate }

      QcConpairAggregate(inputQcConpairAggregate,
                         runQC)
      
      inputFastP4MultiQC
        .join(inputAlfredIgnoreY,by:0)
        .join(inputAlfredIgnoreN,by:0)
        .join(inputConpairConcord4Aggregate,by:0)
        .join(inputConpairContami4Aggregate,by:0)
        .join(inputFacetsQC4CohortMultiQC,by:0)
        .join(inputQualimap4CohortMultiQC,by:0)
        .join(inputHsMetrics, by:0)
        .set{ inputCohortRunMultiQC }

      CohortRunMultiQC(inputCohortRunMultiQC,
                       Channel.value([multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo]),
                       runQC)
    }
  }
}

workflow.onComplete {
  file(params.fileTracking).text = ""
  file(outDir).eachFileRecurse{
    file(params.fileTracking).append(it + "\n")
  }
}