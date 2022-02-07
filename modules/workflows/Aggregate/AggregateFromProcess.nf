include { GermlineAggregateMaf }               from './GermlineAggregateMaf'
include { GermlineAggregateSv }                from './GermlineAggregateSv'
include { QcBamAggregate }                     from './QcBamAggregate'
include { QcConpairAggregate }                 from './QcConpairAggregate'
include { SomaticAggregateFacets }             from './SomaticAggregateFacets'
include { SomaticAggregateLOHHLA }             from './SomaticAggregateLOHHLA'
include { SomaticAggregateMaf }                from './SomaticAggregateMaf'
include { SomaticAggregateMetadata }           from './SomaticAggregateMetadata'
include { SomaticAggregateNetMHC }             from './SomaticAggregateNetMHC'
include { SomaticAggregateSv }                 from './SomaticAggregateSv'
include { CohortRunMultiQC }                   from '../QC/CohortRunMultiQC'

workflow aggregateFromProcess
{
  take:
    inputPairing
    facets4Aggregate
    sv4Aggregate
    snv4Aggregate
    lohhla4Aggregate
    MetaData4Aggregate
    mafFile4AggregateGermline
    sv4AggregateGermline
    bamsQcStats4Aggregate
    collectHsMetricsOutput
    qualimap4Process
    conpair4Aggregate
    FacetsQC4Aggregate
    doWF_facets
    doWF_SV
    doWF_SNV
    doWF_loh
    doWF_mdParse
    doWF_germSV
    doWF_QC
    doWF_germSNV
    fastPJson
    multiqcWesConfig
    multiqcWgsConfig
    multiqcTempoLogo

  main:
    inputPairing.set{ cohortTable }
    cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
        .set{ inputAggregate }

    if (facets4Aggregate){
      input4AggregateFacets = inputAggregate.combine(facets4Aggregate.out.facets4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2],it[4],it[5],it[6],it[7],it[8]]}
    }
    if (sv4Aggregate){
      inputSomaticAggregateSv = inputAggregate.combine(sv4Aggregate.out.sv4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    }
    if (snv4Aggregate){
      inputSomaticAggregateNetMHC = inputAggregate.combine(snv4Aggregate.out.NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputSomaticAggregateMaf    = inputAggregate.combine(snv4Aggregate.out.finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (lohhla4Aggregate){
      inputSomaticAggregateLOHHLA = inputAggregate.combine(lohhla4Aggregate.out.lohhla4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    }
    if (doWF_mdParse){
      inputSomaticAggregateMetadata = inputAggregate.combine(MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_facets && doWF_germSNV){
      inputGermlineAggregateMaf = inputAggregate.combine(mafFile4AggregateGermline, by:[1,2]).groupTuple(by:[2])
    }
    if (sv4AggregateGermline)
    {
      inputGermlineAggregateSv    = inputAggregate.combine(sv4AggregateGermline.out.sv4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique(), it[6].unique()]}
    }

    if (doWF_QC){
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

      if (doWF_QC && params.pairing){
        inputQcConpairAggregate = inputAggregate.combine(conpair4Aggregate.out.conpair4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
        inputAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}.set{ inputFacetsQC4CohortMultiQC }
      }
    }

  if (doWF_SNV){
    SomaticAggregateMaf(inputSomaticAggregateMaf)
    SomaticAggregateFacets(input4AggregateFacets)
    SomaticAggregateNetMHC(inputSomaticAggregateNetMHC)
  }
  if (doWF_SV){
    SomaticAggregateSv(inputSomaticAggregateSv)
  }
  if (doWF_loh){
    SomaticAggregateLOHHLA(inputSomaticAggregateLOHHLA)
  }
  if(doWF_mdParse)
  {
    SomaticAggregateMetadata(inputSomaticAggregateMetadata)
  }
  if (doWF_facets && doWF_germSNV){
    GermlineAggregateMaf(inputGermlineAggregateMaf)
  }
  if (doWF_germSV){
    GermlineAggregateSv(inputGermlineAggregateSv)
  }

  if (doWF_QC && params.pairing){
    inputAlfredIgnoreY.join(inputAlfredIgnoreN)
              .join(inputHsMetrics)
              .set{ inputQcBamAggregate }

    QcBamAggregate(inputQcBamAggregate)
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
