include { SplitLanesR1; SplitLanesR2 } from '../process/Alignment/SplitLanes' 
include { AlignReads }                 from '../process/Alignment/AlignReads'
include { GATK4SPARK_MARKDUPLICATES } from '../nf-core/gatk4spark/markduplicates/main'
include { GATK4_SPLITINTERVALS as BQSR_SPLITINTERVALS } from '../nf-core/gatk4/splitintervals/main'
include { GATK4SPARK_BASERECALIBRATOR } from '../nf-core/gatk4spark/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS } from '../nf-core/gatk4/gatherbqsrreports/main'
include { GATK4SPARK_APPLYBQSR      } from '../nf-core/gatk4spark/applybqsr/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX as INDEX_BQSR_BAM } from '../nf-core/samtools/index/main'

workflow alignment_wf
{
  take:
    inputMapping

  main:

    versions = Channel.empty()

    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    if (params.bamMapping)
    {
      println "Alignment workflow cannot accept bam files for input."
      exit 1
    }
    if(params.mapping)
    {
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
      AlignReads.out.sortedBam
        .groupTuple(by:[3])
        .map { idSample, target, bam, fileID, lane, readIdFile ->
          def idSample_first = idSample instanceof Collection ? idSample.first() : idSample
          def target_first   = target instanceof Collection ? target.first() : target
          if ( !params.watch ){
            for (i in readIdFile.flatten().unique()){
              def readId = "@" + i.getSimpleName().replaceAll("@", ":")
              if(!TempoUtils.checkDuplicates(allReadIds, readId, idSample_first + "\t" + fileID, "the following samples, since they contain the same read: \n${readId}")){exit 1}
            }
          }
          [idSample_first, target_first, bam.flatten().unique()]
        }
        .groupTuple(by: [0])
        .map{ item ->
          def idSample = item[0]
          def target =  item[1] instanceof Collection ? item[1].first() : item[1]
          def bams = item[2].flatten().unique()
          [idSample, bams, target]
        }
        .set { groupedBam }


      GATK4SPARK_MARKDUPLICATES(
	     groupedBam.map{item ->
		def meta = [:]
		meta.id = item[0]
		meta.target = item[2]
		bams = item[1]
		bamSize = 0
		bams.flatten().each{ bamSize = bamSize + it.size()}
		meta.size = bamSize >> 30
		[meta, bams]
	     },
             referenceMap.genomeFile,
             referenceMap.genomeIndex,
	     referenceMap.genomeDict
      )


    // Join with the bai file
      markdup_bam_bai = GATK4_MARKDUPLICATES.out.bam.map{[it[0].id, it[0], it[1]]}
						    .join(GATK4_MARKDUPLICATES.out.bai.map{[it[0].id, it[0], it[1]]}, failOnDuplicate: true, failOnMismatch: true)
						    .map{[it[1], it[2], it[4]]}

      BQSR_SPLITINTERVALS(
	   Channel.from(targetsMap.keySet()).map{ targetId -> [ [ id:"${targetId}"], targetsMap."${targetId}".targetsInterval ]},
	   Channel.fromPath(params.genomes[params.genome].genomeFile).collect().map{ it -> [ [ id:'fasta' ], it ] },
	   Channel.fromPath(params.genomes[params.genome].genomeIndex).collect().map{ it -> [ [ id:'fai' ], it ] },
	   Channel.fromPath(params.genomes[params.genome].genomeDict).collect().map{ it -> [ [ id:'Dict' ], it ] },
      )

      split_interval = BQSR_SPLITINTERVALS.out.split_intervals.map{ item ->
			target = item[0].id
			num_intervals = item[1] instanceof Collection ? item[1].size() : 1
			intervals = item[1]
			[target, num_intervals, intervals]
		}

      bam_and_intervals = markdup_bam_bai.map{[it[0].target, it[0], it[1], it[2]]}
					 .combine(split_interval, by: 0)
					 .transpose()
					 .map{[it[1] + [num_intervals:it[4]], it[2], it[3], it[5]]}
      // Channel Contains [meta, bam, bai, interval]

      GATK4SPARK_BASERECALIBRATOR(
	     bam_and_intervals,
             referenceMap.genomeFile,
             referenceMap.genomeIndex,
             referenceMap.genomeDict,
	     Channel.fromPath(params.genomes[params.genome].knownIndels)
			.concat(Channel.fromPath(params.genomes[params.genome].dbsnp))
			.collect(),
	     Channel.fromPath(params.genomes[params.genome].knownIndelsIndex)
			.concat(Channel.fromPath(params.genomes[params.genome].dbsnpIndex))
			.collect()
      )

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = GATK4SPARK_BASERECALIBRATOR.out.table.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple().branch{
        // Use meta.num_intervals to asses number of intervals
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(table_to_merge.multiple)

    // Mix intervals and no_intervals channels together
    table_bqsr = GATK4_GATHERBQSRREPORTS.out.table.mix(table_to_merge.single.map{ meta, table -> [ meta, table[0] ] })
						  .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }
        // Remove no longer necessary field: num_intervals

    bqsr_input = table_bqsr.map{meta, file -> [meta.id, file]}
			   .combine(bam_and_intervals.map{[it[0].id, it[0], it[1], it[2], it[3]]},by: 0)
			   .map{[it[2], it[3], it[4],  it[1], it[5]]}
    GATK4SPARK_APPLYBQSR(
	bqsr_input,
	referenceMap.genomeFile,
	referenceMap.genomeIndex,
	referenceMap.genomeDict
    )

    bam_to_merge_index = GATK4SPARK_APPLYBQSR.out.bam.map{ meta, bam -> [ groupKey(meta, meta.num_intervals), bam ] }.groupTuple().branch{
        // Use meta.num_intervals to asses number of intervals
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // Only when using intervals
    MERGE_BAM(
	bam_to_merge_index.multiple,
	Channel.fromPath(params.genomes[params.genome].genomeFile).collect().map{ it -> [ [ id:'fasta' ], it ] },
	Channel.fromPath(params.genomes[params.genome].genomeIndex).collect().map{ it -> [ [ id:'fasta_fai' ], it ] }
    )

    // Mix intervals and no_intervals channels together
    bam_all = MERGE_BAM.out.bam.mix(bam_to_merge_index.single.map{ meta, bam -> [ meta, bam[0] ] })
			       .map{ meta, bam -> [ meta - meta.subMap('num_intervals'), bam ] }
        // Remove no longer necessary field: num_intervals

    // Index bam
    INDEX_BQSR_BAM(bam_all)

    // Join with the bai file
    bam_bai = bam_all.join(INDEX_BQSR_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true).map{[it[0].id, it[0].target, it[1], it[2]]}


      File file_bammapping = new File(params.outname)
      file_bammapping.newWriter().withWriter { w ->
          w << "SAMPLE\tTARGET\tBAM\tBAI\n"
      }

      bam_bai
      .map{ idSample, target, bam, bai ->
        [ idSample, target, "${file(params.outDir).toString()}/bams/${idSample}/${idSample}.bam", "${file(params.outDir).toString()}/bams/${idSample}/${idSample}.bam.bai" ]
      }.subscribe { Object obj ->
        file_bammapping.withWriterAppend { out ->
            out.println "${obj[0]}\t${obj[1]}\t${obj[2]}\t${obj[3]}"
        }
      }

    // Gather versions of all tools used
      versions = versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions)
      versions = versions.mix(GATK4SPARK_BASERECALIBRATOR.out.versions)
      versions = versions.mix(GATK4_GATHERBQSRREPORTS.out.versions)
      versions = versions.mix(GATK4SPARK_APPLYBQSR.out.versions)
      versions = versions.mix(MERGE_BAM.out.versions.first())
      versions = versions.mix(INDEX_BQSR_BAM.out.versions.first())
    }
    else{
      if(params.pairing){
        println "ERROR: When --pairing [tsv], --mapping [tsv] must be provided."
        exit 1
      }
    }
  

  emit:
    bam_bai
    fastPJson          = fastPJson
    versions          // channel: [ versions.yml ]
}
