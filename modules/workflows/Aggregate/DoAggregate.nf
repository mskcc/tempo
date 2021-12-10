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

workflow aggregateFromFile
{
  take:
    aggregateFile
    multiqcWesConfig
    multiqcWgsConfig
    multiqcTempoLogo

  main:
    doWF_validate        = true
    doWF_align           = true
    doWF_manta           = true
    doWF_scatter         = true
    doWF_germSNV         = true
    doWF_germSV          = true
    doWF_facets          = true
    doWF_SV              = true
    doWF_loh             = true
    doWF_SNV             = true
    doWF_sampleQC        = true
    doWF_msiSensor       = true
    doWF_mutSig          = true
    doWF_mdParse         = true
    doWF_samplePairingQC = true

    if(!params.watch){
      TempoUtils.extractCohort(file(aggregateFile, checkIfExists: true))
        .set{ inputAggregate }
    }
    else{
      watchAggregateWithPath(file(aggregateFile, checkIfExists: true))
        .set{ inputAggregate }
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

    inputSomaticAggregateMaf      = aggregateList.finalMaf4Aggregate.transpose().groupTuple(by:[2])
    inputSomaticAggregateNetMHC   = aggregateList.NetMhcStats4Aggregate.transpose().groupTuple(by:[2])
    inputPurity4Aggregate         = aggregateList.FacetsPurity4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputHisens4Aggregate         = aggregateList.FacetsHisens4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputOutLog4Aggregate         = aggregateList.FacetsOutLog4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputArmLev4Aggregate         = aggregateList.FacetsArmLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputGeneLev4Aggregate        = aggregateList.FacetsGeneLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputFacetsQC4CohortMultiQC   = aggregateList.FacetsQC4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    inputSomaticAggregateSv       = aggregateList.dellyMantaCombined4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateSvTbi    = aggregateList.dellyMantaCombinedTbi4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputPredictHLA4Aggregate     = aggregateList.predictHLA4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputIntCPN4Aggregate         = aggregateList.intCPN4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateMetadata = aggregateList.MetaData4Aggregate.transpose().groupTuple(by:[2])
    inputGermlineAggregateMaf     = aggregateList.mafFile4AggregateGermline.transpose().groupTuple(by:[2])
    inputGermlineAggregateSv      = aggregateList.dellyMantaCombined4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputGermlineAggregateSvTbi   = aggregateList.dellyMantaCombinedTbi4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputAlfredIgnoreY            = aggregateList.alfredIgnoreYTumor.unique().combine(aggregateList.alfredIgnoreYNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputAlfredIgnoreN            = aggregateList.alfredIgnoreNTumor.unique().combine(aggregateList.alfredIgnoreNNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputQualimap4CohortMultiQC   = aggregateList.qualimapTumor.unique().combine(aggregateList.qualimapNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputHsMetrics                = aggregateList.hsMetricsTumor.unique().combine(aggregateList.hsMetricsNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    aggregateList.conpairConcord4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairConcord4Aggregate}
    aggregateList.conpairContami4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairContami4Aggregate}
    aggregateList.fastpTumor.unique().combine(aggregateList.fastpNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}.set{inputFastP4MultiQC}

    SomaticAggregateMaf(inputSomaticAggregateMaf)

    inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
        .join(inputOutLog4Aggregate, by:[0])
        .join(inputArmLev4Aggregate, by:[0])
        .join(inputGeneLev4Aggregate, by:[0])
        .set{ inputSomaticAggregateFacets }
    SomaticAggregateFacets(inputSomaticAggregateFacets)
    SomaticAggregateNetMHC(inputSomaticAggregateNetMHC)

    inputSomaticAggregateSv.join(inputSomaticAggregateSvTbi)
          .set{ inputSomaticAggregateSv }
    SomaticAggregateSv(inputSomaticAggregateSv)
  
    inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
          .set{ inputSomaticAggregateLOHHLA }
    SomaticAggregateLOHHLA(inputSomaticAggregateLOHHLA)

    SomaticAggregateMetadata(inputSomaticAggregateMetadata)

    GermlineAggregateMaf(inputGermlineAggregateMaf)
  
    inputGermlineAggregateSv.join(inputGermlineAggregateSvTbi)
          .set{ inputGermlineAggregateSv }
    GermlineAggregateSv(inputGermlineAggregateSv)

    inputAlfredIgnoreY.join(inputAlfredIgnoreN)
              .join(inputHsMetrics)
              .set{ inputQcBamAggregate }
    QcBamAggregate(inputQcBamAggregate)

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
      inputSomaticAggregateSv    = inputAggregate.combine(dellyMantaCombined4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputSomaticAggregateSvTbi = inputAggregate.combine(dellyMantaCombinedTbi4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_SNV){
      inputSomaticAggregateNetMHC = inputAggregate.combine(NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputSomaticAggregateMaf    = inputAggregate.combine(finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_loh){
      inputPredictHLA4Aggregate = inputAggregate.combine(predictHLA4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputIntCPN4Aggregate     = inputAggregate.combine(intCPN4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_mdParse){
      inputSomaticAggregateMetadata = inputAggregate.combine(MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])
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
    SomaticAggregateMaf(inputSomaticAggregateMaf)

    inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
        .join(inputOutLog4Aggregate, by:[0])
        .join(inputArmLev4Aggregate, by:[0])
        .join(inputGeneLev4Aggregate, by:[0])
        .set{ inputSomaticAggregateFacets }

    SomaticAggregateFacets(inputSomaticAggregateFacets)
    SomaticAggregateNetMHC(inputSomaticAggregateNetMHC)
  }
  if (doWF_SV){
    inputSomaticAggregateSv.join(inputSomaticAggregateSvTbi)
          .set{ inputSomaticAggregateSv }

    SomaticAggregateSv(inputSomaticAggregateSv)
  }
  if (doWF_loh){
    inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
          .set{ inputSomaticAggregateLOHHLA }

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
