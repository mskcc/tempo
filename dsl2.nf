#!/usr/bin/env nextflow
nextflow.enable.dsl=2


if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println "ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)"
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
wallTimeExitCode = params.wallTimeExitCode ? params.wallTimeExitCode.split(',').collect{it.trim().toLowerCase()} : []
multiqcWesConfig = workflow.projectDir + "/lib/multiqc_config/exome_multiqc_config.yaml"
multiqcWgsConfig = workflow.projectDir + "/lib/multiqc_config/wgs_multiqc_config.yaml"
multiqcTempoLogo = workflow.projectDir + "/docs/tempoLogo.png"
epochMap = [:]
startEpoch = new Date().getTime()
limitInputLines = 0
chunkSizeLimit = params.chunkSizeLimit
if (params.watch == true){
  touchInputs()
}

include { defineReferenceMap; loadTargetReferences  } from './modules/local/define_maps'

include { SplitLanesR1; SplitLanesR2 } from './modules/process/SplitLanes' addParams(wallTimeExitCode: wallTimeExitCode)
include { AlignReads }                 from './modules/process/AlignReads' addParams(outDir: outDir, wallTimeExitCode: wallTimeExitCode)


println ""

pairingQc = params.pairing

if (!runSomatic && runGermline){
    println "WARNING: You can't run GERMLINE section without running SOMATIC section. Activating SOMATIC section automatically"
    runSomatic = true
}

if (runAggregate == false){
  if (!params.mapping && !params.bamMapping){
    println "ERROR: (--mapping/-bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv] ) or (--aggregate [tsv]) need to be provided, otherwise nothing to be run."
    exit 1
  }
}
else if (runAggregate == true){
  if ((params.mapping || params.bamMapping) && params.pairing){
    if (!(runSomatic || runGermline || runQC)){
      println "ERROR: Nothing to be aggregated. One or more of the option --somatic/--germline/--QC need to be enabled when using --aggregate"
    }
  }
  else if ((params.mapping || params.bamMapping) && !params.pairing){
    if (!runQC){
      println "ERROR: Nothing to be aggregated. --QC need to be enabled when using --mapping/--bamMapping [tsv], --pairing false and --aggregate true."
      exit 1
    }
  }
  else{
    println "ERROR: (--mapping/--bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv]) or (--aggregate [tsv]) need to be provided when using --aggregate true"
    println "       If you want to run aggregate only, you need to use --aggregate [tsv]. See manual"
    exit 1
  }

}
else {
  if ((runSomatic || runGermline || runQC) && !params.mapping && !params.bamMapping){
    println "ERROR: Conflict input! When running --aggregate [tsv] with --mapping/--bamMapping/--pairing [tsv] disabled, --QC/--somatic/--germline all need to be disabled!"
    println "       If you want to run aggregate somatic/germline/qc, just include an additianl colum PATH in the [tsv] and no need to use --QC/--somatic/--germline flag, since it's auto detected. See manual"
    exit 1
  }
}

if (!(params.cosmic in ['v2', 'v3'])) {
  println "ERROR: Possible values of mutational signature reference --cosmic is 'v2', 'v3'"
  exit 1
}


workflow {
  referenceMap = defineReferenceMap()
  targetsMap   = loadTargetReferences()

  mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)
  inputMapping = TempoUtils.extractFastq(mappingFile, params.assayType, targetsMap.keySet())

  // Skip these processes if starting from aligned BAM files
  if (params.mapping) {

    // Parse input FASTQ mapping
    if(params.watch != true){
    inputMapping.groupTuple(by: [0])
                .map { idSample, targets, files_pe1, files_pe2
                  -> tuple(groupKey(idSample, targets.size()), targets, files_pe1, files_pe2)
                }
                .transpose()
          .set{ inputMapping }
    }

    inputMapping.map{ idSample, target, file_pe1, file_pe2 ->
                    [idSample, target, file_pe1, file_pe2, idSample + "@" + file_pe1.getSimpleName(), file_pe1.getSimpleName()]
                }
                .set{ inputFastqs }

    if (params.splitLanes) {
        inputFastqs.set{ fastqsNeedSplit }
        inputFastqs.set{ fastqsNoNeedSplit }

        fastqsNeedSplit
            .filter{ item -> !(item[2].getName() =~ /_L(\d){3}_/) }
            .multiMap{ idSample, target, file_pe1, file_pe2, fileID, lane ->
                inputFastqR1: [idSample, target, file_pe1, file_pe1.toString()]
                inputFastqR2: [idSample, target, file_pe2, file_pe2.toString()]
            }
            .set{ fastqsNeedSplit }

        fastqsNoNeedSplit
            .filter{ item -> item[2].getName() =~ /_L(\d){3}_/ }
            .map { idSample, target, file_pe1, file_pe2, fileID, lane
                    -> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
            }
            .set{ fastqsNoNeedSplit }

        perLaneFastqsR1 = SplitLanesR1(fastqsNeedSplit.inputFastqR1).R1SplitData
        perLaneFastqsR2 = SplitLanesR2(fastqsNeedSplit.inputFastqR2).R2SplitData


        def fastqR1fileIDs = [:]
        perLaneFastqsR1 = perLaneFastqsR1.transpose()
            .map{ item ->
                def idSample = item[0]
                def target = item[1]
                def fastq = item[2]
                def fileID = idSample + "@" + item[3].getSimpleName()
                def lane = fastq.getSimpleName().split("_L00")[1].split("_")[0]
                def laneCount = item[4].getSimpleName().toInteger()

                // This only checks if same read groups appears in two or more fastq files which belongs to the same sample. Cross sample check will be performed after AlignReads since the read group info is not available for fastqs which does not need to be split.
                if ( !params.watch ){
                if(!TempoUtils.checkDuplicates(fastqR1fileIDs, fileID + "@" + lane, fileID + "\t" + fastq, "the follwoing fastq files since they contain the same RGID")){exit 1}
                }

                [idSample, target, fastq, fileID, lane, laneCount]
            }

        def fastqR2fileIDs = [:]
        perLaneFastqsR2 = perLaneFastqsR2.transpose()
                .map{ item ->
                def idSample = item[0]
                def target = item[1]
                def fastq = item[2]
                def fileID = idSample + "@" + item[3].getSimpleName()
                def lane = fastq.getSimpleName().split("_L00")[1].split("_")[0]
                def laneCount = item[4].getSimpleName().toInteger()

                if ( !params.watch ){
                    if(!TempoUtils.checkDuplicates(fastqR2fileIDs, fileID + "@" + lane, fileID + "\t" + fastq, "the follwoing fastq files since they contain the same RGID")){exit 1}
                }

                [idSample, target, fastq, fileID, lane, laneCount]
                }

        fastqFiles  = perLaneFastqsR1
            .mix(perLaneFastqsR2)
            .groupTuple(by: [0,1,3,4,5], size: 2, sort: true)
            .map {  idSample, target, fastqPairs, fileID, lanes, laneCount ->
            tuple(idSample, target, fastqPairs, groupKey(fileID, laneCount), lanes)
            }
            .map{ idSample, target, fastqPairs, fileID, lane ->
                [idSample, target, fastqPairs[0], fastqPairs[1], fileID, lane]
            }
            .map{ item ->
                    def idSample = item[0]
                    def target = item[1]
                    def fastqPair1 = item[2]
                    def fastqPair2 = item[3]
                    if (item[2].toString().split("_R1").size() < item[3].toString().split("_R1").size()) {
                    fastqPair1 = item[3]
                    fastqPair2 = item[2]
                    }
                    def fileID = item[4]
                    def lane = item[5]

                [idSample, target, fastqPair1, fastqPair1.size(), fastqPair2, fastqPair2.size(), fileID, lane]
            }
            .mix(fastqsNoNeedSplit)
    }
    else{
     fastqFiles = inputFastqs.map { idSample, target, file_pe1, file_pe2, fileID, lane
					-> tuple(idSample, target, file_pe1, file_pe1.size(), file_pe2, file_pe2.size(), groupKey(fileID, 1), lane)
				   }
    }


    //Align reads to reference.
    readAlignOut = AlignReads(fastqFiles, Channel.value([referenceMap.genomeFile, referenceMap.bwaIndex]))

    fastPJson4cohortMultiQC = readAlignOut.fastPJson4MultiQC
    .groupTuple(by:[2])
    .map{idSample, jsonFile, fileID -> 
      def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
      [idSampleout, jsonFile]
    }.groupTuple(by: [0])
    .map{ idSample, jsonFile -> 
      [idSample, jsonFile.flatten()]
    }
  
    fastPJson4sampleMultiQC = readAlignOut.fastPJson4MultiQC
    .groupTuple(by:[2])
    .map{idSample, jsonFile, fileID -> 
      def idSampleout = idSample[0] instanceof Collection ? idSample[0].first() : idSample[0]
      [idSampleout, jsonFile]
    }.groupTuple(by: [0])
    .map{ idSample, jsonFile -> 
      [idSample, jsonFile.flatten()]
    }

    // Check for FASTQ files which might have different path but contains the same reads, based only on the name of the first read.
    def allReadIds = [:]
    readAlignOut.sortedBam.map { idSample, target, bam, fileID, lane, readIdFile -> def readId = "@" + readIdFile.getSimpleName().replaceAll("@", ":")
      // Use the first line of the fastq file (the name of the first read) as unique identifier to check across all the samples if there is any two fastq files contains the same read name, if so, we consider there are some human error of mixing up the same reads into different fastq files
      if ( !params.watch ){
       if(!TempoUtils.checkDuplicates(allReadIds, readId, idSample + "\t" + bam, "the following samples, since they contain the same read: \n${readId}")){exit 1}
      }
      [idSample, target, bam, fileID, lane]
	   }
	   .groupTuple(by: [3])
	   .map{ item ->
		      def idSample = item[0] instanceof Collection ? item[0].first() : item[0]
		      def target   = item[1] instanceof Collection ? item[1].first() : item[1]
		      def bams = item[2]
		      [idSample, target, bams]
	   }
	   .groupTuple(by: [0])
	   .map{ item ->
		      def idSample = item[0]
		      def target =  item[1] instanceof Collection ? item[1].first() : item[1]
		      def bams = item[2].flatten()
		      [idSample, bams, target]
	   }
	   .set{ groupedBam }

  } //end if params.mapping

}