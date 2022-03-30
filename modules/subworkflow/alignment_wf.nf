include { SplitLanesR1; SplitLanesR2 } from '../process/Alignment/SplitLanes' 
include { AlignReads }                 from '../process/Alignment/AlignReads'
include { MergeBamsAndMarkDuplicates } from '../process/Alignment/MergeBamsAndMarkDuplicates'
include { RunBQSR }                    from '../process/Alignment/RunBQSR' 

workflow alignment_wf
{
  take:
    inputMapping

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

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
            [idSample, target, file_pe1, file_pe2, idSample + '@' + file_pe1.getSimpleName(), file_pe2.getSimpleName()]
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


      File file = new File(params.outname)
      file.newWriter().withWriter { w ->
          w << "SAMPLE\tTARGET\tBAM\tBAI\n"
      }

      RunBQSR.out.bamResults.subscribe { Object obj ->
        file.withWriterAppend { out ->
            out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}"
        }
      }
  

  emit:
    RunBQSR_bamsBQSR   = RunBQSR.out.bamsBQSR
    RunBQSR_bamResults = RunBQSR.out.bamResults
    RunBQSR_bamSize    = RunBQSR.out.bamSize
    fastPJson          = fastPJson
}
