include { GermlineAggregateMaf }               from './GermlineAggregateMaf'
include { GermlineAggregateSv }                from './GermlineAggregateSv'
include { QcBamAggregate }                     from './QcBamAggregate'
include { QcConpairAggregate }                 from './QcConpairAggregate'
include { AggregateFacets }                    from './AggregateFacets'
include { AggregateLOHHLA }                    from './AggregateLOHHLA'
include { AggregateMaf }                       from './AggregateMaf'
include { AggregateMetadata }                  from './AggregateMetadata'
include { AggregateNetMHC }                    from './AggregateNetMHC'
include { AggregateSv }                        from './AggregateSv'
include { CohortRunMultiQC }                   from '../QC/CohortRunMultiQC'

workflow aggregateFromProcess
{
  take:
    inputPairing
    FacetsPurity4Aggregate
    FacetsHisens4Aggregate
    FacetsOutLog4Aggregate
    FacetsArmLev4Aggregate
    FacetsGeneLev4Aggregate
    dellyMantaCombined4Aggregate
    dellyMantaCombinedTbi4Aggregate
    NetMhcStats4Aggregate
    finalMaf4Aggregate
    predictHLA4Aggregate
    intCPN4Aggregate
    MetaData4Aggregate
    mafFile4AggregateGermline
    dellyMantaCombined4AggregateGermline
    dellyMantaCombinedTbi4AggregateGermline
    bamsQcStats4Aggregate
    collectHsMetricsOutput
    qualimap4Process
    conpairConcord4Aggregate
    conpairContami4Aggregate
    FacetsQC4Aggregate
    doWF_facets
    doWF_SV
    doWF_SNV
    doWF_loh
    doWF_mdParse
    doWF_germSV
    doWF_sampleQC
    doWF_germSNV
    doWF_samplePairingQC
    fastPJson
    multiqcWesConfig
    multiqcWgsConfig
    multiqcTempoLogo

  main:
    inputPairing.set{ cohortTable }
    cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
        .set{ inputAggregate }

    if (doWF_facets){
      inputPurity4Aggregate    = inputAggregate.combine(FacetsPurity4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputHisens4Aggregate    = inputAggregate.combine(FacetsHisens4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputOutLog4Aggregate    = inputAggregate.combine(FacetsOutLog4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputArmLev4Aggregate    = inputAggregate.combine(FacetsArmLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputGeneLev4Aggregate   = inputAggregate.combine(FacetsGeneLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_SV){
      inputAggregateSv    = inputAggregate.combine(dellyMantaCombined4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputAggregateSvTbi = inputAggregate.combine(dellyMantaCombinedTbi4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_SNV){
      inputAggregateNetMHC = inputAggregate.combine(NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputAggregateMaf    = inputAggregate.combine(finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_loh){
      inputPredictHLA4Aggregate = inputAggregate.combine(predictHLA4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputIntCPN4Aggregate     = inputAggregate.combine(intCPN4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_mdParse){
      inputAggregateMetadata = inputAggregate.combine(MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_facets && doWF_germSNV){
      inputGermlineAggregateMaf = inputAggregate.combine(mafFile4AggregateGermline, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_germSV)
    {
      inputGermlineAggregateSv    = inputAggregate.combine(dellyMantaCombined4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
      inputGermlineAggregateSvTbi = inputAggregate.combine(dellyMantaCombinedTbi4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    }

    if (doWF_sampleQC){
      bamsQcStats4Aggregate.branch{ item ->
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
      
      inputPairing.combine(collectHsMetricsOutput)
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

      inputPairing.combine(qualimap4Process)
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

      if (doWF_samplePairingQC){
        inputAggregate.combine(conpairConcord4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairConcord4Aggregate }
        inputAggregate.combine(conpairContami4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairContami4Aggregate }
        inputAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}.set{ inputFacetsQC4CohortMultiQC }
      }
    }

  if (doWF_SNV){
    AggregateMaf(inputAggregateMaf)

    inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
        .join(inputOutLog4Aggregate, by:[0])
        .join(inputArmLev4Aggregate, by:[0])
        .join(inputGeneLev4Aggregate, by:[0])
        .set{ inputAggregateFacets }

    AggregateFacets(inputAggregateFacets)
    AggregateNetMHC(inputAggregateNetMHC)
  }
  if (doWF_SV){
    inputAggregateSv.join(inputAggregateSvTbi)
          .set{ inputAggregateSv }

    AggregateSv(inputAggregateSv)
  }
  if (doWF_loh){
    inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
          .set{ inputAggregateLOHHLA }

    AggregateLOHHLA(inputAggregateLOHHLA)
  }
  if(doWF_mdParse)
  {
    AggregateMetadata(inputAggregateMetadata)
  }
  if (doWF_facets && doWF_germSNV){
    GermlineAggregateMaf(inputGermlineAggregateMaf)
  }
  if (doWF_germSV){
    // --- Aggregate per-sample germline data, SVs
    inputGermlineAggregateSv.join(inputGermlineAggregateSvTbi)
          .set{ inputGermlineAggregateSv }

    GermlineAggregateSv(inputGermlineAggregateSv)
  }

  if (doWF_sampleQC){
    inputAlfredIgnoreY.join(inputAlfredIgnoreN)
              .join(inputHsMetrics)
              .set{ inputQcBamAggregate }

    QcBamAggregate(inputQcBamAggregate)

    if(doWF_samplePairingQC)
    {
      inputConpairConcord4Aggregate.join(inputConpairContami4Aggregate)
          .set{ inputQcConpairAggregate }

      QcConpairAggregate(inputQcConpairAggregate)
      
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
                        Channel.value([multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo]))
    }
  }
}
