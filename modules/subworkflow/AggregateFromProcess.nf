include { GermlineAggregateMaf }               from '../process/Aggregate/GermlineAggregateMaf'
include { GermlineAggregateSv }                from '../process/Aggregate/GermlineAggregateSv'
include { QcBamAggregate }                     from '../process/Aggregate/QcBamAggregate'
include { QcConpairAggregate }                 from '../process/Aggregate/QcConpairAggregate'
include { SomaticAggregateFacets }             from '../process/Aggregate/SomaticAggregateFacets'
include { SomaticAggregateLOHHLA }             from '../process/Aggregate/SomaticAggregateLOHHLA'
include { SomaticAggregateMaf }                from '../process/Aggregate/SomaticAggregateMaf'
include { SomaticAggregateMetadata }           from '../process/Aggregate/SomaticAggregateMetadata'
include { SomaticAggregateNetMHC }             from '../process/Aggregate/SomaticAggregateNetMHC'
include { SomaticAggregateSv }                 from '../process/Aggregate/SomaticAggregateSv'
include { SomaticAggregateSvSignatures }       from '../process/Aggregate/SomaticAggregateSvSignatures'
include { SomaticAggregateHRDetect }           from '../process/Aggregate/SomaticAggregateHRDetect'
include { SomaticAggregateSVclone }            from '../process/Aggregate/SomaticAggregateSVclone'
include { CohortRunMultiQC }                   from '../process/Aggregate/CohortRunMultiQC'
include { watchMapping; watchBamMapping; watchPairing; watchAggregateWithResult; watchAggregate } from '../function/watch_inputs.nf'

workflow aggregateFromProcess
{
  take:
    inputPairing
    runAggregate
    facets4Aggregate
    sv4Aggregate
    snv4Aggregate
    hrd4Aggregate
    svclone4Aggregate
    lohhla4Aggregate
    MetaData4Aggregate
    snv4AggregateGermline
    sv4AggregateGermline
    sampleQC4Aggregate
    conpair4Aggregate
    fastPJson
    multiqcWesConfig
    multiqcWgsConfig
    multiqcTempoLogo

  main:
    if (runAggregate != true){
      if (!params.watch){
        TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
	          .groupTuple()
		  .map{ cohort, idTumor, idNormal, pathNoUse
		        -> tuple( groupKey(cohort, idTumor instanceof Collection ? idTumor.size() : 1), idTumor, idNormal)
		  }
		  .transpose()
		  .set{inputAggregate}
      }
      else{
        watchAggregate(file(runAggregate, checkIfExists: false))
	              .set{inputAggregate}
      }
    }
    else {
      inputPairing.set{ cohortTable }
      cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
                 .set{ inputAggregate }
    }

    if (facets4Aggregate){
      input4AggregateFacets = inputAggregate.combine(facets4Aggregate.out.facets4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2],it[4],it[5],it[6],it[7],it[8]]}
    }
    if (sv4Aggregate){
      inputSomaticAggregateSv = 
        inputAggregate.combine(sv4Aggregate.out.sv4Aggregate, by:[1,2])
          .groupTuple(by:[2])
          .map{[it[2], it[4]]}

      inputSomaticAggregateSvSignatures = 
        inputAggregate.combine(sv4Aggregate.out.SVSignatures, by:[1,2])
          .groupTuple(by:[2])
          .map{[it[2], it[4], it[5]]}
    }
    if (snv4Aggregate){
      inputSomaticAggregateNetMHC = inputAggregate.combine(snv4Aggregate.out.NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputSomaticAggregateMaf    = inputAggregate.combine(snv4Aggregate.out.finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (hrd4Aggregate){
      inputSomaticAggregateHrd = inputAggregate.combine(
        hrd4Aggregate.out.HRDetect, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (svclone4Aggregate){
      inputSomaticAggregateSVclone = inputAggregate
        .combine(svclone4Aggregate.out.svclone4Aggregate, by:[1,2])
        .groupTuple(by:[2])
        .map{[it[2], it[4], it[5]]}
    }
    if (lohhla4Aggregate){
      inputSomaticAggregateLOHHLA = inputAggregate.combine(lohhla4Aggregate.out.lohhla4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    }
    if (MetaData4Aggregate){
      inputSomaticAggregateMetadata = inputAggregate.combine(MetaData4Aggregate.out.MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (snv4AggregateGermline){
      inputGermlineAggregateMaf = inputAggregate.combine(snv4AggregateGermline.out.snv4AggregateGermline, by:[1,2]).groupTuple(by:[2])
    }
    if (sv4AggregateGermline)
    {
      inputGermlineAggregateSv = inputAggregate.combine(sv4AggregateGermline.out.sv4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    }

    if (sampleQC4Aggregate){
      sampleQC4Aggregate.out.bamsQcStats4Aggregate.branch{ item ->
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
      
      inputPairing.combine(sampleQC4Aggregate.out.collectHsMetricsOutput)
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

      inputPairing.combine(sampleQC4Aggregate.out.qualimap4Process)
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
    }

    if (conpair4Aggregate){
      inputQcConpairAggregate = inputAggregate.combine(conpair4Aggregate.out.conpair4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
      FacetsQC4Aggregate = facets4Aggregate ? facets4Aggregate.out.FacetsQC4Aggregate : inputPairing.map{ idTumor, idNormal -> ["placeHolder",idTumor, idNormal,"",""]}
      inputFacetsQC4CohortMultiQC = inputAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    }

  if (facets4Aggregate){
    SomaticAggregateFacets(input4AggregateFacets)
  }
  if (snv4Aggregate){
    SomaticAggregateMaf(inputSomaticAggregateMaf)
    SomaticAggregateNetMHC(inputSomaticAggregateNetMHC)
  }
  if (sv4Aggregate){
    SomaticAggregateSv(inputSomaticAggregateSv)
    if (params.assayType == "genome"){
      SomaticAggregateSvSignatures(inputSomaticAggregateSvSignatures)
    }
  }
  if (hrd4Aggregate){
    SomaticAggregateHRDetect(inputSomaticAggregateHrd)
  }
  if (svclone4Aggregate){
    SomaticAggregateSVclone(inputSomaticAggregateSVclone)
  }
  if (lohhla4Aggregate){
    SomaticAggregateLOHHLA(inputSomaticAggregateLOHHLA)
  }
  if(MetaData4Aggregate)
  {
    SomaticAggregateMetadata(inputSomaticAggregateMetadata)
  }
  if (snv4AggregateGermline){
    GermlineAggregateMaf(inputGermlineAggregateMaf)
  }
  if (sv4AggregateGermline){
    GermlineAggregateSv(inputGermlineAggregateSv)
  }

  if (sampleQC4Aggregate){
    inputAlfredIgnoreY.join(inputAlfredIgnoreN)
              .join(inputHsMetrics)
              .set{ inputQcBamAggregate }

    QcBamAggregate(inputQcBamAggregate)
  }

  if (conpair4Aggregate) {
    QcConpairAggregate(inputQcConpairAggregate)

    inputFastP4MultiQC
      .join(inputAlfredIgnoreY,by:0)
      .join(inputAlfredIgnoreN,by:0)
      .join(inputQcConpairAggregate,by:0)
      .join(inputFacetsQC4CohortMultiQC,by:0)
      .join(inputQualimap4CohortMultiQC,by:0)
      .join(inputHsMetrics, by:0)
      .set{ inputCohortRunMultiQC }

    CohortRunMultiQC(inputCohortRunMultiQC,
                      Channel.value([multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo]))
  }
}