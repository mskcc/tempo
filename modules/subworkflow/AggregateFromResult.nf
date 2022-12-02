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
include { CohortRunMultiQC }                   from '../process/Aggregate/CohortRunMultiQC'
include { watchMapping; watchBamMapping; watchPairing; watchAggregateWithResult; watchAggregate } from '../function/watch_inputs.nf'

workflow aggregateFromResult
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
      watchAggregateWithResult(file(aggregateFile, checkIfExists: true))
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
      sv4Aggregate: params.assayType == "genome" ?
        [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.combined.annot.filtered.pass.clustered.bedpe")] :
        [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.combined.annot.filtered.bedpe")]
      svSignatures4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*_exposures.tsv"), file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*_catalogues.pdf")]
      hrd4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.hrdetect.tsv")]
      predictHLA4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.DNA.HLAlossPrediction_CI.txt")]
      intCPN4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*DNA.IntegerCPN_CI.txt")]
      MetaData4Aggregate: [idTumor, idNormal, cohort, "placeHolder", file(path + "/somatic/" + idTumor + "__" + idNormal + "/*/*.sample_data.txt")]
      mafFile4AggregateGermline: [idTumor, idNormal, cohort, "placeHolder", file(path + "/germline/" + idNormal + "/*/" + idTumor + "__" + idNormal + ".germline.final.maf")]
      sv4AggregateGermline: [idNormal, cohort, idTumor, "placeHolder", "noTumor", file(path + "/germline/" + idNormal + "/*/*.combined.annot.filtered.pass.bedpe")]
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
    inputSomaticAggregateSv       = aggregateList.sv4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateSvSignatures = aggregateList.svSignatures4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    inputSomaticAggregateHrd      = aggregateList.hrd4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputPredictHLA4Aggregate     = aggregateList.predictHLA4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputIntCPN4Aggregate         = aggregateList.intCPN4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateMetadata = aggregateList.MetaData4Aggregate.transpose().groupTuple(by:[2])
    inputGermlineAggregateMaf     = aggregateList.mafFile4AggregateGermline.transpose().groupTuple(by:[2])
    inputGermlineAggregateSv      = aggregateList.sv4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
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

    SomaticAggregateSv(inputSomaticAggregateSv)
    SomaticAggregateSvSignatures(inputSomaticAggregateSvSignatures)
    SomaticAggregateHRDetect(inputSomaticAggregateHrd)
  
    inputPredictHLA4Aggregate.join(inputIntCPN4Aggregate)
          .set{ inputSomaticAggregateLOHHLA }
    SomaticAggregateLOHHLA(inputSomaticAggregateLOHHLA)

    SomaticAggregateMetadata(inputSomaticAggregateMetadata)

    GermlineAggregateMaf(inputGermlineAggregateMaf)
  
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
