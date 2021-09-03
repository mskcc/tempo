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
include { SplitLanesR1; SplitLanesR2 } from './modules/process/Alignment/SplitLanes' addParams(wallTimeExitCode: wallTimeExitCode)
include { AlignReads }                 from './modules/process/Alignment/AlignReads' addParams(outDir: outDir, wallTimeExitCode: wallTimeExitCode)
include { MergeBamsAndMarkDuplicates } from './modules/process/Alignment/MergeBamsAndMarkDuplicates' addParams(wallTimeExitCode: wallTimeExitCode)
include { RunBQSR }                    from './modules/process/Alignment/RunBQSR' addParams(outDir: outDir, wallTimeExitCode: wallTimeExitCode)
include { CrossValidateSamples }       from './modules/process/QC/SampleValidation'

//Somatic
include { CreateScatteredIntervals }   from './modules/process/Somatic/CreateScatteredIntervals'
include { RunMutect2 }                 from './modules/process/Somatic/RunMutect2'
include { SomaticCombineMutect2Vcf }   from './modules/process/Somatic/SomaticCombineMutect2Vcf' addParams(outDir: outDir)
include { SomaticDellyCall }           from './modules/process/Somatic/SomaticDellyCall' addParams(outDir: outDir)
include { SomaticRunManta }            from './modules/process/Somatic/SomaticRunManta' addParams(outDir: outDir)
include { SomaticMergeDellyAndManta }  from './modules/process/Somatic/SomaticMergeDellyAndManta' addParams(outDir: outDir)
include { SomaticRunStrelka2 }         from './modules/process/Somatic/SomaticRunStrelka2' addParams(outDir: outDir)
include { SomaticCombineChannel }      from './modules/process/Somatic/SomaticCombineChannel' addParams(outDir: outDir)
include { SomaticAnnotateMaf }         from './modules/process/Somatic/SomaticAnnotateMaf' addParams(outDir: outDir)
include { RunMutationSignatures }      from './modules/process/Somatic/RunMutationSignatures'
include { DoFacets }                   from './modules/process/Somatic/DoFacets' addParams(outDir: outDir)
include { DoFacetsPreviewQC }          from './modules/process/Somatic/DoFacetsPreviewQC' addParams(outDir: outDir)
include { RunPolysolver }              from './modules/process/Somatic/RunPolysolver'
include { RunLOHHLA }                  from './modules/process/Somatic/RunLOHHLA' addParams(outDir: outDir)
include { RunNeoantigen }              from './modules/process/Somatic/RunNeoantigen' addParams(outDir: outDir, wallTimeExitCode: wallTimeExitCode)
include { SomaticFacetsAnnotation }    from './modules/process/Somatic/SomaticFacetsAnnotation' addParams(outDir: outDir)
include { RunMsiSensor }               from './modules/process/Somatic/RunMsiSensor'
include { MetaDataParser }             from './modules/process/Somatic/MetaDataParser' addParams(outDir: outDir)

//Germline
include { GermlineRunHaplotypecaller }        from './modules/process/Germline/GermlineRunHaplotypecaller'
include { GermlineCombineHaplotypecallerVcf } from './modules/process/Germline/GermlineCombineHaplotypecallerVcf' addParams(outDir: outDir)
include { GermlineRunStrelka2 }               from './modules/process/Germline/GermlineRunStrelka2' addParams(outDir: outDir)
include { GermlineCombineChannel }            from './modules/process/Germline/GermlineCombineChannel' addParams(outDir: outDir)
include { GermlineAnnotateMaf }               from './modules/process/Germline/GermlineAnnotateMaf' addParams(outDir: outDir)
include { GermlineFacetsAnnotation }          from './modules/process/Germline/GermlineFacetsAnnotation' addParams(outDir: outDir)
include { GermlineDellyCall }                 from './modules/process/Germline/GermlineDellyCall' addParams(outDir: outDir)
include { GermlineRunManta }                  from './modules/process/Germline/GermlineRunManta' addParams(outDir: outDir)
include { GermlineMergeDellyAndManta }        from './modules/process/Germline/GermlineMergeDellyAndManta' addParams(outDir: outDir)

//QC
include { QcCollectHsMetrics }                 from './modules/process/QC/QcCollectHsMetrics' addParams(outDir: outDir, wallTimeExitCode: wallTimeExitCode)
include { QcQualimap }                         from './modules/process/QC/QcQualimap' addParams(outDir: outDir)
include { QcAlfred }                           from './modules/process/QC/QcAlfred' addParams(outDir: outDir)
include { SampleRunMultiQC }                   from './modules/process/QC/SampleRunMultiQC'
include { QcPileup }                           from './modules/process/QC/QcPileup'
include { QcConpair }                          from './modules/process/QC/QcConpair'
include { SomaticRunMultiQC }                  from './modules/process/QC/SomaticRunMultiQC'
include { QcConpairAll }                       from './modules/process/QC/QcConpairAll'

println ''

pairingQc = params.pairing


workflow {
  referenceMap = defineReferenceMap()
  targetsMap   = loadTargetReferences()

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
        //(checkPairing1, checkPairing2, inputPairing) = TempoUtils.extractPairing(pairingFile).into(3)
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
      if (!runSomatic && !runGermline && !runQC){
        println "ERROR: --pairing [tsv] is not used because none of --somatic/--germline/--QC was enabled. If you only need to do BAM QC and/or BAM generation, remove --pairing [tsv]."
        exit 1
      }
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
  if (params.mapping) {
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

    fastPJson4cohortMultiQC = AlignReads.out.fastPJson4MultiQC
      .groupTuple(by:[2])
      .map { idSample, jsonFile, fileID ->
        def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
        [idSampleout, jsonFile]
      }.groupTuple(by: [0])
      .map { idSample, jsonFile ->
        [idSample, jsonFile.flatten()]
      }

    fastPJson4sampleMultiQC = AlignReads.out.fastPJson4MultiQC
      .groupTuple(by:[2])
      .map { idSample, jsonFile, fileID ->
        def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
        [idSampleout, jsonFile]
      }.groupTuple(by: [0])
      .map { idSample, jsonFile ->
        [idSample, jsonFile.flatten()]
      }

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
  } //end if params.mapping


  /*
  ================================================================================
  =                                PAIRING TUMOR and NORMAL                      =
  ================================================================================
  */
  // If starting with BAM files, parse BAM pairing input
  if (params.bamMapping) {
    //inputMapping.into{bamsBQSR4Alfred; bamsBQSR4Qualimap; bamsBQSR4CollectHsMetrics; bamsBQSR4Tumor; bamsBQSR4Normal; bamsBQSR4QcPileup; bamPaths4MultiQC}
    if (runQC){
      inputMapping.map{idSample, target, bam, bai ->
        [ idSample,target, bam.getParent() ]
      }.set{locateFastP4MultiQC}

      inputMapping.map{ idSample,target, bamFolder -> 
        [idSample, file(bamFolder + "/fastp/*json")]
      }.set{fastPJson}
    } 
  }

  if (params.pairing) {
    // Parse input FASTQ mapping
    //inputPairing.into{pairing4T; pairing4N; pairingTN; pairing4QC; inputPairing}
    RunBQSR.out.bamsBQSR.combine(inputPairing)
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

    RunBQSR.out.bamsBQSR.combine(inputPairing)
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
        [ targetId, targetsMap[targetId]."targetsBedGz", targetsMap[targetId]."targetsBedGzTbi" ]
      }

  targets4Intervals = Channel.from(targetsMap.keySet())
      .map{ targetId ->
        [ targetId, targetsMap."${targetId}"."targetsBedGz", targetsMap."${targetId}"."targetsBedGzTbi" ]
      }

    CreateScatteredIntervals(Channel.value([referenceMap.genomeFile, 
                                            referenceMap.genomeIndex, 
                                            referenceMap.genomeDict]), 
                                            targets4Intervals,
                                            runSomatic, runGermline)

    //Associating interval_list files with BAM files, putting them into one channel
    //bamFiles.into{bamsTN4Intervals; bamsForDelly; bamsForManta; bams4Strelka; bamns4CombineChannel; bamsForMsiSensor; bamFiles4DoFacets; bamsForLOHHLA }
    bamFiles.combine(CreateScatteredIntervals.out.mergedIList, by: 2).map{
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
    RunBQSR.out.bamsBQSR.map{ idSample, target, bam, bai ->
        [idSample, target, bam, bai, targetsMap."$target".targetsInterval,  targetsMap."$target".baitsInterval]
    }.set{bamsBQSR4HsMetrics}

    QcCollectHsMetrics(bamsBQSR4HsMetrics,
                       Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
                       params.assayType,
                       runQC)

    if (runQC && params.assayType != "exome"){
      QcCollectHsMetrics.out.collectHsMetricsOutput
        .map{ idSample, target, bam, bai, targetList, baitList -> [idSample, ""]}
        .into{ collectHsMetricsOutput }
    }

    RunBQSR.out.bamsBQSR
      .map{ idSample, target, bam, bai -> [ idSample, target, bam, bai, file(targetsMap."$target".targetsBed) ]}
      .set{ bamsBQSR4Qualimap }

    QcQualimap(bamsBQSR4Qualimap,
               runQC)

    Channel.from(true, false).set{ ignore_read_groups }
    RunBQSR.out.bamsBQSR
      .map{ idSample, target, bam, bai -> 
        [ idSample, target, bam, bai, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi ]
      }.set{ bamsBQSR4Alfred }

    QcAlfred(ignore_read_groups, 
             bamsBQSR4Alfred,
             Channel.value([referenceMap.genomeFile]),
             runQC)

    QcAlfred.out.alfredOutput
      .groupTuple(size:2, by:0)
      .join(fastPJson4sampleMultiQC, by:0)
      .join(QcQualimap.out.qualimap4Process, by:0)
      .join(QcCollectHsMetrics.out.collectHsMetricsOutput, by:0)
      .set{ sampleMetrics4MultiQC }


    SampleRunMultiQC(sampleMetrics4MultiQC, 
                     Channel.value([multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo]))


    if (pairingQc) {
      QcPileup(RunBQSR.out.bamsBQSR,
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

}
