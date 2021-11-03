#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println 'ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)'
  exit 1
}

// User-set runtime parameters
outDir           = file(params.outDir).toAbsolutePath()
outname          = params.outname
runGermline      = params.germline
runSomatic       = params.somatic
runQC            = params.QC
runAggregate     = params.aggregate
runConpairAll    = false
wallTimeExitCode = params.wallTimeExitCode ? params.wallTimeExitCode.split(',').collect { it.trim().toLowerCase() } : []
multiqcWesConfig = workflow.projectDir + '/lib/multiqc_config/exome_multiqc_config.yaml'
multiqcWgsConfig = workflow.projectDir + '/lib/multiqc_config/wgs_multiqc_config.yaml'
multiqcTempoLogo = workflow.projectDir + '/docs/tempoLogo.png'
limitInputLines  = 0
chunkSizeLimit   = params.chunkSizeLimit

//Utility Includes
include { defineReferenceMap; loadTargetReferences } from './modules/local/define_maps'
include { touchInputs; watchMapping; watchBamMapping; watchPairing; watchAggregateWithPath; watchAggregate } from './modules/local/watch_inputs'

//Aggregate Includes
include { GermlineAggregateMaf }               from './modules/workflows/Aggregate/GermlineAggregateMaf'
include { GermlineAggregateSv }                from './modules/workflows/Aggregate/GermlineAggregateSv'
include { QcBamAggregate }                     from './modules/workflows/Aggregate/QcBamAggregate'
include { QcConpairAggregate }                 from './modules/workflows/Aggregate/QcConpairAggregate'
include { SomaticAggregateFacets }             from './modules/workflows/Aggregate/SomaticAggregateFacets'
include { SomaticAggregateLOHHLA }             from './modules/workflows/Aggregate/SomaticAggregateLOHHLA'
include { SomaticAggregateMaf }                from './modules/workflows/Aggregate/SomaticAggregateMaf'
include { SomaticAggregateMetadata }           from './modules/workflows/Aggregate/SomaticAggregateMetadata'
include { SomaticAggregateNetMHC }             from './modules/workflows/Aggregate/SomaticAggregateNetMHC'
include { SomaticAggregateSv }                 from './modules/workflows/Aggregate/SomaticAggregateSv'

pairingQc    = params.pairing
referenceMap = defineReferenceMap()
targetsMap   = loadTargetReferences()

//Sub-workflow Includes
include { validate_wf }        from './modules/workflows/WorkflowControls/validate_wf'         addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { alignment_wf }       from './modules/workflows/WorkflowControls/alignment_wf'        addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { manta_wf }           from './modules/workflows/WorkflowControls/manta_wf'            addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { msiSensor_wf }       from './modules/workflows/WorkflowControls/msiSensor_wf'        addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { mutSig_wf }          from './modules/workflows/WorkflowControls/mutSig_wf' 
include { mdParse_wf }         from './modules/workflows/WorkflowControls/mdParse_wf' 
include { loh_wf }             from './modules/workflows/WorkflowControls/loh_wf'              addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { facets_wf }          from './modules/workflows/WorkflowControls/facets_wf'           addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { sv_wf }              from './modules/workflows/WorkflowControls/sv_wf'               addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { snv_wf }             from './modules/workflows/WorkflowControls/snv_wf'              addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { sampleQC_wf }        from './modules/workflows/WorkflowControls/sampleQC_wf'         addParams(referenceMap: referenceMap, targetsMap: targetsMap, multiqcWesConfig: multiqcWesConfig, multiqcWgsConfig: multiqcWgsConfig, multiqcTempoLogo: multiqcTempoLogo)
include { samplePairingQC_wf } from './modules/workflows/WorkflowControls/samplePairingQC_wf'  addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { somaticMultiQC_wf }  from './modules/workflows/WorkflowControls/somaticMultiQC_wf'   addParams(multiqcWesConfig: multiqcWesConfig, multiqcWgsConfig: multiqcWgsConfig, multiqcTempoLogo: multiqcTempoLogo)
include { scatter_wf }         from './modules/workflows/WorkflowControls/scatter_wf'          addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { germlineSNV_wf }     from './modules/workflows/WorkflowControls/germlineSNV_wf'      addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { germlineSNV_facets } from './modules/workflows/WorkflowControls/germlineSNV_facets' 
include { germlineSV_wf }      from './modules/workflows/WorkflowControls/germlineSV_wf'       addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { PairTumorNormal }    from './modules/workflows/WorkflowControls/PairTumorNormal' 


workflow {
  //Set flags for when each pipeline is required to run.
  doWF_validate        = (params.pairing || params.bamMapping) ? true : false
  doWF_align           = (params.mapping) ? true : false
  doWF_manta           = (params.snvWF || params.svWF || params.mutsigWF ) ? true : false
  doWF_scatter         = (params.snvWF || params.mutsigWF || params.germSNV) ? true : false
  doWF_germSNV         = (params.germSNV) ? true : false
  doWF_germSV          = (params.germSV) ? true : false
  doWF_facets          = (params.lohWF || params.facetsWF || params.snvWF || params.mutsigWF || params.germSNV) ? true : false
  doWF_SV              = (params.svWF) ? true : false
  doWF_loh             = (params.lohWF || params.snvWF || params.mutsigWF) ? true : false
  doWF_SNV             = (params.snvWF || params.mutsigWF) ? true : false
  doWF_sampleQC        = (params.sampleQCWF || params.samplePairingQCWF) ? true : false
  doWF_msiSensor       = (params.msiWF) ? true : false
  doWF_mutSig          = (params.mutsigWF) ? true : false
  doWF_mdParse         = (doWF_manta && doWF_scatter && doWF_facets && doWF_loh && doWF_SNV && doWF_msiSensor && doWF_mutSig) ? true : false
  doWF_samplePairingQC = (params.samplePairingQCWF) ? true : false

  //Handle specific invalid conditions.
  if (runAggregate == false) {
    if (!params.mapping && !params.bamMapping) {
      println 'ERROR: (--mapping/-bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv] ) or (--aggregate [tsv]) need to be provided, otherwise nothing to be run.'
      exit 1
    }
  }
  else if (runAggregate == true) {
    if ((params.mapping || params.bamMapping) && params.pairing) {
      if (!(runSomatic || runGermline || runQC)) {
        println 'ERROR: Nothing to be aggregated. One or more of the option --somatic/--germline/--QC need to be enabled when using --aggregate'
      }
    }
    else if ((params.mapping || params.bamMapping) && !params.pairing) {
      if (!params.sampleQCWF || params.samplePairingQCWF) {
        println 'ERROR: Nothing to be aggregated. --sampleQCWF or --samplePairingQCWF need to be enabled when using --mapping/--bamMapping [tsv], --pairing false and --aggregate true.'
        exit 1
      }
    }
    else {
      println 'ERROR: (--mapping/--bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv]) or (--aggregate [tsv]) need to be provided when using --aggregate true'
      println '       If you want to run aggregate only, you need to use --aggregate [tsv]. See manual'
      exit 1
    }
  }
  else {
    if ((runSomatic || runGermline || runQC) && !params.mapping && !params.bamMapping) {
      println 'ERROR: Conflict input! When running --aggregate [tsv] with --mapping/--bamMapping/--pairing [tsv] disabled, --QC/--somatic/--germline all need to be disabled!'
      println "       If you want to run aggregate somatic/germline/qc, just include an additianl colum PATH in the [tsv] and no need to use --QC/--somatic/--germline flag, since it's auto detected. See manual"
      exit 1
    }
  }

  //Align can run without pairing/cross validate, but only if nothing else is running.
  if((!params.pairing) && (doWF_validate || doWF_manta || doWF_scatter || doWF_germSNV || doWF_germSV 
                    || doWF_facets || doWF_SV || doWF_loh || doWF_SNV || doWF_sampleQC 
                    || doWF_msiSensor || doWF_mutSig || doWF_mdParse || doWF_samplePairingQC))
  {
    println "ERROR: Certain workflows cannot be performed without pairing information."
    println "\tProvide a --pairing file, or disable other sub-workflows to proceed."
    exit 1
  }
  else
  {
    if (params.watch == false) {
      keySet = targetsMap.keySet()
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)      
      inputMapping = params.mapping ? TempoUtils.extractFastq(mappingFile, params.assayType, targetsMap.keySet()) : TempoUtils.extractBAM(mappingFile, params.assayType, targetsMap.keySet())
    }
    else if (params.watch == true) {
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: false) : file(params.bamMapping, checkIfExists: false)
      inputMapping  = params.mapping ? watchMapping(mappingFile, params.assayType, targetsMap.keySet()) : watchBamMapping(mappingFile, params.assayType, targetsMap.keySet())
    }
  }

  //Begin executing modules for the run.
  if (doWF_validate) {
    validate_wf()
    inputMapping = validate_wf.out.inputMapping
    inputPairing = validate_wf.out.inputPairing
  }
  else{
    if(params.pairing){
      println "ERROR: When --pairing [tsv], --mapping/--bamMapping [tsv] must be provided."
      exit 1
    }
  }

  if (doWF_align)
  {
    alignment_wf(inputMapping)
  }

  //Handle input bams as coming originally from bams, or from an alignment this run.
  if (params.bamMapping) {
    inputBam = inputMapping
    if (doWF_sampleQC || doWF_samplePairingQC){
      inputMapping.map{idSample, target, bam, bai ->
        [ idSample,target, bam.getParent() ]
      }.set{ locateFastP4MultiQC }

      locateFastP4MultiQC.map{ idSample,target, bamFolder -> 
        [idSample, file(bamFolder + "/fastp/*json")]
      }.set{ fastPJson }
    } 
  }
  else
  {
    inputBam = alignment_wf.out.RunBQSR_bamsBQSR
    fastPJson = alignment_wf.out.fastPJson
  }

  //Generate several channels for tumor/normal pairing required for downstream processes.
  if (params.pairing) {
    PairTumorNormal(inputBam, inputPairing)
    bamFiles   = PairTumorNormal.out.bamFiles
    bams       = PairTumorNormal.out.bams
    bamsNormal = PairTumorNormal.out.bamsNormal
    bamsTumor  = PairTumorNormal.out.bamsTumor
  } 

  if(doWF_manta)
  {
    manta_wf(bamFiles)
  }

  if(doWF_scatter)
  {
    scatter_wf()
  }

  if(doWF_germSNV)
  {
    germlineSNV_wf(bams, bamsTumor, scatter_wf.out.mergedIList)
  }

  if(doWF_germSV)
  {
    germlineSV_wf(bams)
  }

  if(doWF_facets)
  {
    facets_wf(bamFiles)
    if(doWF_germSNV)
    {
      germlineSNV_facets(facets_wf.out.facetsForMafAnno, germlineSNV_wf.out.mafFileGermline)
    }
  }

  if(doWF_SV)
  {
    sv_wf(bamFiles, manta_wf.out.manta4Combine)
  }

  if(doWF_loh)
  {
    loh_wf(bams, bamFiles, facets_wf.out.facetsPurity)
  }

  if(doWF_SNV)
  {
    snv_wf(bamFiles, scatter_wf.out.mergedIList, manta_wf.out.mantaToStrelka, loh_wf.out.hlaOutput, facets_wf.out.facetsForMafAnno)
  }

  if(doWF_sampleQC)
  {
    sampleQC_wf(inputBam, fastPJson)
  }

  if(doWF_msiSensor)
  {
    msiSensor_wf(bamFiles)
  }

  if(doWF_mutSig)
  {
    mutSig_wf(snv_wf.out.mafFile)
  }

  if(doWF_mdParse)
  {
    facets_wf.out.facetsPurity.combine(snv_wf.out.maf4MetaDataParser, by: [0,1,2])
      .combine(facets_wf.out.FacetsQC4MetaDataParser, by: [0,1,2])
      .combine(msiSensor_wf.out.msi4MetaDataParser, by: [0,1,2])
      .combine(mutSig_wf.out.mutSig4MetaDataParser, by: [0,1,2])
      .combine(loh_wf.out.hlaOutput, by: [1,2])
      .unique()
      .map{ idNormal, target, idTumor, purityOut, mafFile, qcOutput, msifile, mutSig, placeHolder, polysolverFile ->
      [idNormal, target, idTumor, purityOut, mafFile, qcOutput, msifile, mutSig, placeHolder, polysolverFile, targetsMap."$target".codingBed]
    }.set{ mergedChannelMetaDataParser }

    mdParse_wf(mergedChannelMetaDataParser)
  }

  if(doWF_samplePairingQC)
  {
    samplePairingQC_wf(inputBam, inputPairing, runConpairAll)

    if (doWF_SNV) {
      samplePairingQC_wf.out.conpairOutput
        .map{ placeHolder, idTumor, idNormal, conpairFiles -> [idTumor, idNormal, conpairFiles]}
        .join(facets_wf.out.FacetsQC4SomaticMultiQC, by:[0,1])
        .set{ somaticMultiQCinput }
      
      FacetsQC4Aggregate = facets_wf.out.FacetsQC4Aggregate
    } else {
      samplePairingQC_wf.out.conpairOutput
        .map{ placeHolder, idTumor, idNormal, conpairFiles -> [idTumor, idNormal, conpairFiles, "", ""]}
        .set{ somaticMultiQCinput }

      inputPairing.map{ idTumor, idNormal -> 
        ["placeHolder",idTumor, idNormal,"",""]
      }.set{ FacetsQC4Aggregate }
    }
    somaticMultiQC_wf(somaticMultiQCinput)
  }
  //This needs to be done so that 'FacetsQC4Aggregate' exists during aggregation for any circumstance.
  else
  {
    if(!doWF_validate)
    {
      inputPairing = Channel.empty()
      inputPairing.map{ idTumor, idNormal -> 
        ["placeHolder","", "","",""]
      }.set{ FacetsQC4Aggregate }
    }
    else{
      inputPairing.map{ idTumor, idNormal -> 
        ["placeHolder",idTumor, idNormal,"",""]
      }.set{ FacetsQC4Aggregate }
    }
  }
  


/*
================================================================================
=                              Cohort Aggregation                              =
================================================================================
*/
  if ( !params.mapping && !params.bamMapping ){
    runSomatic = true
    runGermline = true
    runQC = true
    pairingQc = true
    if(!params.watch){
      TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
        .set{ inputAggregate }
    }
    else{
      watchAggregateWithPath(file(runAggregate, checkIfExists: true))
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

    inputSomaticAggregateMaf = aggregateList.finalMaf4Aggregate.transpose().groupTuple(by:[2])
    inputSomaticAggregateNetMHC = aggregateList.NetMhcStats4Aggregate.transpose().groupTuple(by:[2])
    inputPurity4Aggregate = aggregateList.FacetsPurity4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputHisens4Aggregate = aggregateList.FacetsHisens4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputOutLog4Aggregate = aggregateList.FacetsOutLog4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    finalMaf4Aggregate = aggregateList.FacetsArmLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputGeneLev4Aggregate = aggregateList.FacetsGeneLev4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputFacetsQC4CohortMultiQC = aggregateList.FacetsQC4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4], it[5]]}
    inputSomaticAggregateSv = aggregateList.dellyMantaCombined4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputPredictHLA4Aggregate = aggregateList.predictHLA4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputIntCPN4Aggregate = aggregateList.intCPN4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}
    inputSomaticAggregateMetadata = aggregateList.MetaData4Aggregate.transpose().groupTuple(by:[2])
    inputGermlineAggregateMaf = aggregateList.mafFile4AggregateGermline.transpose().groupTuple(by:[2])
    inputGermlineAggregateSv = aggregateList.dellyMantaCombined4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputGermlineAggregateSvTbi = aggregateList.dellyMantaCombinedTbi4AggregateGermline.transpose().groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    inputAlfredIgnoreY = aggregateList.alfredIgnoreYTumor.unique().combine(aggregateList.alfredIgnoreYNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputAlfredIgnoreN = aggregateList.alfredIgnoreNTumor.unique().combine(aggregateList.alfredIgnoreNNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputQualimap4CohortMultiQC = aggregateList.qualimapTumor.unique().combine(aggregateList.qualimapNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    inputHsMetrics = aggregateList.hsMetricsTumor.unique().combine(aggregateList.hsMetricsNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}
    aggregateList.conpairConcord4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairConcord}
    aggregateList.conpairContami4Aggregate.transpose().groupTuple(by:[2]).map{[it[2], it[4]]}.set{inputConpairContami}
    aggregateList.fastpTumor.unique().combine(aggregateList.fastpNormal.unique(), by:[0,1,2]).transpose().groupTuple(by:[0]).map{ [it[0], it[3].unique(), it[4].unique()]}.set{inputFastP4MultiQC}
  
  }
  else if(!(runAggregate == false)) {
    if (!(runAggregate == true)){
      if(!params.watch){
        TempoUtils.extractCohort(file(runAggregate, checkIfExists: true))
            .groupTuple()
            .map{ cohort, idTumor, idNormal, pathNoUse
                  -> tuple( groupKey(cohort, idTumor instanceof Collection ? idTumor.size() : 1), idTumor, idNormal)
            }
            .transpose()
            .set{ inputAggregate }
      }
      else{
        watchAggregate(file(runAggregate, checkIfExists: false))
          .set{ inputAggregate }
      }
    }
    else if(runAggregate == true){
      inputPairing.set{ cohortTable }
      cohortTable.map{ idTumor, idNormal -> ["default_cohort", idTumor, idNormal]}
          .set{ inputAggregate }
    }
    else{}

    if (doWF_facets){
      inputPurity4Aggregate    = inputAggregate.combine(facets_wf.out.FacetsPurity4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputHisens4Aggregate    = inputAggregate.combine(facets_wf.out.FacetsHisens4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputOutLog4Aggregate    = inputAggregate.combine(facets_wf.out.FacetsOutLog4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputArmLev4Aggregate    = inputAggregate.combine(facets_wf.out.FacetsArmLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputGeneLev4Aggregate   = inputAggregate.combine(facets_wf.out.FacetsGeneLev4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_SV){
      inputSomaticAggregateSv    = inputAggregate.combine(sv_wf.out.dellyMantaCombined4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputSomaticAggregateSvTbi = inputAggregate.combine(sv_wf.out.dellyMantaCombinedTbi4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_SNV){
      inputSomaticAggregateNetMHC = inputAggregate.combine(snv_wf.out.NetMhcStats4Aggregate, by:[1,2]).groupTuple(by:[2])
      inputSomaticAggregateMaf    = inputAggregate.combine(snv_wf.out.finalMaf4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_loh){
      inputPredictHLA4Aggregate = inputAggregate.combine(loh_wf.out.predictHLA4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
      inputIntCPN4Aggregate     = inputAggregate.combine(loh_wf.out.intCPN4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}
    }
    if (doWF_mdParse){
      inputSomaticAggregateMetadata = inputAggregate.combine(mdParse_wf.out.MetaData4Aggregate, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_facets && doWF_germSNV){
      inputGermlineAggregateMaf = inputAggregate.combine(germlineSNV_facets.out.mafFile4AggregateGermline, by:[1,2]).groupTuple(by:[2])
    }
    if (doWF_germSV)
    {
      inputGermlineAggregateSv    = inputAggregate.combine(germlineSV_wf.out.dellyMantaCombined4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
      inputGermlineAggregateSvTbi = inputAggregate.combine(germlineSV_wf.out.dellyMantaCombinedTbi4AggregateGermline, by:[2]).groupTuple(by:[1]).map{[it[1], it[5].unique()]}
    }

    if (doWF_sampleQC){
      sampleQC_wf.out.bamsQcStats4Aggregate.branch{ item ->
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
      
      inputPairing.combine(sampleQC_wf.out.collectHsMetricsOutput)
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

      inputPairing.combine(sampleQC_wf.out.qualimap4Process)
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
        inputAggregate.combine(samplePairingQC_wf.out.conpairConcord4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairConcord4Aggregate }
        inputAggregate.combine(samplePairingQC_wf.out.conpairContami4Aggregate, by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4]]}.set{ inputConpairContami4Aggregate }
        inputAggregate.combine(FacetsQC4Aggregate,by:[1,2]).groupTuple(by:[2]).map{[it[2], it[4], it[5]]}.set{ inputFacetsQC4CohortMultiQC }
      }
    }
  }
  else{}

  if (params.watch == true) {
    epochMap = [:]
    for (i in ["mapping","bamMapping","pairing","aggregate"]) {
      if (params.containsKey(i)){ epochMap[params."${i}"] = 0 }
    }
    startEpoch = new Date().getTime()
    touchInputs(chunkSizeLimit, startEpoch, epochMap)
  }

  if (runAggregate) {
    if (doWF_facets){
      SomaticAggregateMaf(inputSomaticAggregateMaf)

      inputPurity4Aggregate.join(inputHisens4Aggregate, by:[0])
		     .join(inputOutLog4Aggregate, by:[0])
		     .join(inputArmLev4Aggregate, by:[0])
		     .join(inputGeneLev4Aggregate, by:[0])
		     .set{ inputSomaticAggregateFacets }

      SomaticAggregateFacets(inputSomaticAggregateFacets)
    }
    if (doWF_SNV){
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
}
workflow.onComplete {
  file(params.fileTracking).text = ""
  file(outDir).eachFileRecurse{
    file(params.fileTracking).append(it + "\n")
  }
}
