def touchInputs() {
  new Timer().schedule({
  def timeNow = new Date().getTime()
  limitInputLines = chunkSizeLimit + ( ((timeNow - startEpoch)/60000) * (chunkSizeLimit / params.touchInputsInterval) )
  for ( i in epochMap.keySet() ){
    fileEpoch = file(i).lastModified()
    if (( fileEpoch > epochMap[i]) || (chunkSizeLimit > 0 )) {
      epochMap[i] = fileEpoch
      "touch -ca ${i}".execute()
    }
  }
} as TimerTask, 1000, params.touchInputsInterval * 60 * 1000 ) // convert minutes to milliseconds
}

def watchMapping(tsvFile, assayType) {
  def index = 1 
  Channel.watchPath( tsvFile, 'create, modify' )
	 .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
	 .view()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def fastqFile1 = file(row.FASTQ_PE1, checkIfExists: false)
              def fastqFile2 = file(row.FASTQ_PE2, checkIfExists: false)
              def numOfPairs = row.NUM_OF_PAIRS.toInteger()
              if(!TempoUtils.checkTarget(target, assayType, targetsMap.keySet())){}
              if(!TempoUtils.checkNumberOfItem(row, 5, tsvFile)){}

              [idSample, numOfPairs, target, fastqFile1, fastqFile2]
      }
      .map{ idSample, numOfPairs, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, numOfPairs), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchBamMapping(tsvFile, assayType){
  def index = 1
  Channel.watchPath( tsvFile, 'create, modify' )
	 .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
	 .map{ row ->
              def idSample = row.SAMPLE
              def target = row.TARGET
              def bam = file(row.BAM, checkIfExists: false)
              def bai = file(row.BAI, checkIfExists: false)
              if(!TempoUtils.checkTarget(target, assayType, targetsMap.keySet())){}
              if(!TempoUtils.checkNumberOfItem(row, 4, tsvFile)){}

              [idSample, target, bam, bai]
      }
      .map{ idSample, target, files_pe1, files_pe2
              -> tuple( groupKey(idSample, 1), target, files_pe1, files_pe2)
      }
      .transpose()
      .unique()
}

def watchPairing(tsvFile){
  Channel.watchPath( tsvFile, 'create, modify' )
	 .splitCsv(sep: '\t', header: true)
	 .unique()
         .map { row ->
              def TUMOR_ID = row.TUMOR_ID
              def NORMAL_ID = row.NORMAL_ID
              if(!TempoUtils.checkNumberOfItem(row, 2, tsvFile)){}

              [TUMOR_ID, NORMAL_ID]
         }
	 .unique()
}

def watchAggregateWithPath(tsvFile) {
  def index = 1 
  Channel.watchPath(tsvFile, 'create, modify')
     .map{ row -> 
	      index = 1 
	      row
	 }.splitCsv(sep: '\t', header: true)
	 .map{ row -> 
	      [index++] + row
	 }.filter{ row ->
	      if (chunkSizeLimit > 0 ){
	      	row[0] < limitInputLines
	      } else { 1 }
	 }.map{ row ->
	      row[1]
	 }.unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              def path = row.PATH
              if(!TempoUtils.checkNumberOfItem(row, 5, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal, path]
         }
         .map { cohort, cohortSize, idTumor, idNormal, path
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal, path)
         }
         .transpose()
	 .unique()
}

def watchAggregate(tsvFile) {
  Channel.watchPath(file(runAggregate), 'create, modify')
     .splitCsv(sep: '\t', header: true)
	 .unique()
         .map{ row ->
              def idNormal = row.NORMAL_ID
              def idTumor = row.TUMOR_ID
              def cohort = row.COHORT
              def cohortSize = row.COHORT_SIZE.toInteger()
              if(!TempoUtils.checkNumberOfItem(row, 4, file(runAggregate))){}

              [cohort, cohortSize, idTumor, idNormal]
         }
         .map { cohort, cohortSize, idTumor, idNormal
                    -> tuple( groupKey(cohort, cohortSize), idTumor, idNormal)
         }
         .transpose()
	 .unique()
}

