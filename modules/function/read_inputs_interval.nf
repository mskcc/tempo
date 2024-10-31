workflow watchMapping {
  take:
    read_inputs_channel
  main:
    def index = 0
    def interval_count = 0
    read_inputs_channel
	.map{
          interval_count = it + 1
          index = 0
          file(params.mapping)
        }
	.splitCsv(sep: '\t', header: true)
	.filter{ row ->
	    index += 1
	    if (params.chunkSizeLimit > 0 ){
	        index <= params.chunkSizeLimit*interval_count
	    }else{
	        1
	    }
	}.unique()
	.map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: false)
              def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: false)
              def numOfPairs = row.NUM_OF_PAIRS.toInteger()
              if(!TempoUtils.checkTarget(target, params.assayType, validTargetsList)){}
              if(!TempoUtils.checkNumberOfItem(row, 5, params.mapping)){}

              [idSample, numOfPairs, target, fastqFile1, fastqFile2]
         }
         .map{ idSample, numOfPairs, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, numOfPairs), target, files_pe1, files_pe2)
         }
         .transpose()
         .unique()
	 .set{mapping_ch}
  emit:
    mapping_ch
}

workflow watchBamMapping {
  take:
  read_inputs_channel

  main:
  def index = 0
  def interval_count = 0
  read_inputs_channel
      .map{ 
          interval_count = it + 1
	  index = 0
          file(params.bamMapping)
      }
      .splitCsv(sep: '\t', header: true)
      .filter{ row ->
          index = index + 1
          if (params.chunkSizeLimit > 0 ){
              index <= params.chunkSizeLimit*interval_count
          }else{ 1 }
	  
      }.unique()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def bam = file(row.BAM, checkIfExists: false)
              def bai = file(row.BAI, checkIfExists: false)
              if(!TempoUtils.checkTarget(target, params.assayType, params.targetsMap.keySet())){}
              if(!TempoUtils.checkNumberOfItem(row, 4, params.bamMapping)){}

              [idSample, target, bam, bai]
      }
      .map{ idSample, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, 1), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
      .set{bamMapping_ch}
   emit:
     bamMapping_ch

}

workflow watchPairing {
  take:
    read_inputs_channel
  main:
    read_inputs_channel
      .map{ file(params.pairing) }
      .splitCsv(sep: '\t', header: true)
      .unique()
      .map { row ->
         def TUMOR_ID = row.TUMOR_ID
         def NORMAL_ID = row.NORMAL_ID
         if(!TempoUtils.checkNumberOfItem(row, 2, params.pairing)){}

         [TUMOR_ID, NORMAL_ID]
      }.unique()
      .set{pairing_ch}
  emit:
    pairing_ch

}

workflow watchAggregateWithResult {
  take:
    read_inputs_channel
  main:
    def index = 0
    def interval_count = 0

    read_inputs_channel
      .map{
        interval_count = it + 1
	index = 0
        file(params.aggregate)
      }.splitCsv(sep: '\t', header: true)
      .filter{ row ->
        index += 1
        if (params.chunkSizeLimit > 0 ){
	  index <= params.chunkSizeLimit*interval_count
	} else { 1 }
    }.map{ row ->
      def idNormal = row.NORMAL_ID
      def idTumor = row.TUMOR_ID
      def cohort = row.COHORT
      def cohortSize = row.COHORT_SIZE.toInteger()
      def path = row.PATH
      if(!TempoUtils.checkNumberOfItem(row, 5, file(params.aggregate))){}
      [cohort, cohortSize, idTumor, idNormal, path]
    }.map { cohort, cohortSize, idTumor, idNormal, path
      -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal, path)
    }.transpose()
    .unique()
    .set{aggregate_ch}

  emit:
    aggregate_ch

}

workflow watchAggregate {
  take:
    read_inputs_channel
  main:
    read_inputs_channel
        .map{ file(params.aggregate) }
     .splitCsv(sep: '\t', header: true)
	 .unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              if(!TempoUtils.checkNumberOfItem(row, 4, file(params.aggregate))){}

              [cohort, cohortSize, idTumor, idNormal]
         }
         .map { cohort, cohortSize, idTumor, idNormal
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal)
         }
         .transpose()
	 .unique()
	 .set{aggregate_ch}
  emit:
  aggregate_ch
}
