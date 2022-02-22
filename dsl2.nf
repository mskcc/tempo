#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

if (!(workflow.profile in ['juno', 'awsbatch', 'docker', 'singularity', 'test_singularity', 'test'])) {
  println 'ERROR: You need to set -profile (values: juno, awsbatch, docker, singularity)'
  exit 1
}

// User-set runtime parameters
outDir           = file(params.outDir).toAbsolutePath()
outname          = params.outname
runAggregate     = params.aggregate
runConpairAll    = false
wallTimeExitCode = params.wallTimeExitCode ? params.wallTimeExitCode.split(',').collect { it.trim().toLowerCase() } : []
multiqcWesConfig = workflow.projectDir + '/lib/multiqc_config/exome_multiqc_config.yaml'
multiqcWgsConfig = workflow.projectDir + '/lib/multiqc_config/wgs_multiqc_config.yaml'
multiqcTempoLogo = workflow.projectDir + '/docs/tempoLogo.png'
params.startEpoch = new Date().getTime() 


//Utility Includes
include { defineReferenceMap; loadTargetReferences } from './modules/local/define_maps'
include { touchInputs; watchMapping; watchBamMapping; watchPairing; watchAggregateWithResult; watchAggregate } from './modules/local/watch_inputs'

pairingQc    = params.pairing
referenceMap = defineReferenceMap()
targetsMap   = loadTargetReferences()

//Sub-workflow Includes
include { validate_wf }          from './modules/workflows/WorkflowControls/validate_wf'         addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { alignment_wf }         from './modules/workflows/WorkflowControls/alignment_wf'        addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { manta_wf }             from './modules/workflows/WorkflowControls/manta_wf'            addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { msiSensor_wf }         from './modules/workflows/WorkflowControls/msiSensor_wf'        addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { mutSig_wf }            from './modules/workflows/WorkflowControls/mutSig_wf' 
include { mdParse_wf }           from './modules/workflows/WorkflowControls/mdParse_wf' 
include { loh_wf }               from './modules/workflows/WorkflowControls/loh_wf'              addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { facets_wf }            from './modules/workflows/WorkflowControls/facets_wf'           addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { sv_wf }                from './modules/workflows/WorkflowControls/sv_wf'               addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { snv_wf }               from './modules/workflows/WorkflowControls/snv_wf'              addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { sampleQC_wf }          from './modules/workflows/WorkflowControls/sampleQC_wf'         addParams(referenceMap: referenceMap, targetsMap: targetsMap, multiqcWesConfig: multiqcWesConfig, multiqcWgsConfig: multiqcWgsConfig, multiqcTempoLogo: multiqcTempoLogo)
include { samplePairingQC_wf }   from './modules/workflows/WorkflowControls/samplePairingQC_wf'  addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { somaticMultiQC_wf }    from './modules/workflows/WorkflowControls/somaticMultiQC_wf'   addParams(multiqcWesConfig: multiqcWesConfig, multiqcWgsConfig: multiqcWgsConfig, multiqcTempoLogo: multiqcTempoLogo)
include { scatter_wf }           from './modules/workflows/WorkflowControls/scatter_wf'          addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { germlineSNV_wf }       from './modules/workflows/WorkflowControls/germlineSNV_wf'      addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { germlineSV_wf }        from './modules/workflows/WorkflowControls/germlineSV_wf'       addParams(referenceMap: referenceMap, targetsMap: targetsMap)
include { PairTumorNormal }      from './modules/workflows/WorkflowControls/PairTumorNormal' 
include { aggregateFromResult }    from './modules/workflows/Aggregate/AggregateFromResult'
include { aggregateFromProcess } from './modules/workflows/Aggregate/AggregateFromProcess'
include { hrdetect_wf }          from './modules/workflows/WorkflowControls/hrdetect_wf'
include { ascat_wf }             from './modules/workflows/WorkflowControls/ascat_wf'            addParams(referenceMap: referenceMap)

aggregateParamIsFile = !(runAggregate instanceof Boolean)
// check if --aggregate is a file

WFs = params.workflows instanceof Boolean ? '' : params.workflows

WFs = WFs.split(',').collect{it.trim().toLowerCase()}.unique()

WFs = (!params.mapping && !params.bamMapping && aggregateParamIsFile) ? ['snv','sv','mutsig','germSNV','germSV','lohhla','facets','qc','msisensor'] : WFs

workflow {
  //Set flags for when each pipeline is required to run.
  doWF_align           = (params.mapping) ? true : false
  doWF_manta           = ['snv', 'sv', 'mutsig'].any(it -> it in WFs) ? true : false
  doWF_scatter         = ['snv', 'sv', 'mutsig', 'germsnv'].any(it -> it in WFs) ? true : false
  doWF_germSNV         = 'germsnv' in WFs ? true : false
  doWF_germSV          = 'germsv' in WFs ? true : false
  doWF_facets          = ['lohhla', 'facets', 'snv', 'mutsig', 'germsnv'].any(it -> it in WFs) ? true : false
  doWF_SV              = 'sv' in WFs ? true : false
  doWF_loh             = ['lohhla', 'snv', 'mutsig'].any(it -> it in WFs) ? true : false
  doWF_SNV             = ['snv', 'mutsig'].any(it -> it in WFs) ? true : false ? true : false
  doWF_QC	       = 'qc' in WFs ? true : false
  doWF_msiSensor       = 'msisensor' in WFs ? true : false
  doWF_mutSig          = 'mutsig' in WFs ? true : false
  doWF_mdParse         = (doWF_manta && doWF_scatter && doWF_facets && doWF_loh && doWF_SNV && doWF_msiSensor && doWF_mutSig) ? true : false

  doWF_AggregateFromResult = false
  doWF_AggregateFromProcess = false

  if (!params.mapping && !params.bamMapping) {
      if (aggregateParamIsFile) { doWF_AggregateFromResult = true }
      else {
        println 'ERROR: (--mapping/-bamMapping [tsv]) or (--mapping/--bamMapping [tsv] & --pairing [tsv] ) or (--aggregate [tsv]) need to be provided, otherwise nothing to be run.'
        exit 1
      }
  }
  else if (WFs == [''] && runAggregate) {
        println 'ERROR: No provided sub-workflows enabled to aggregate. Remove --aggregate or add sub-workflows"'
        exit 1
  }
  else{
      doWF_AggregateFromProcess = runAggregate ? true : false
  }

  if (!params.pairing && WFs != ['qc'] && WFs != [''] && !doWF_AggregateFromResult){
      println "ERROR: Certain workflows cannot be performed without pairing information."
      println "\tProvide a --pairing [tsv], or disable other sub-workflows to proceed."
      exit 1
  }

  if (params.bamMapping && WFs == ['']){
      println "ERROR: No sub-workflows to run.."
      println "\tPlease provide sub-workflows using --workflow parameters."
      exit 1
  }

  if(params.pairing && !params.mapping && !params.bamMapping){
    println "ERROR: When --pairing [tsv], --mapping/--bamMapping [tsv] must be provided."
    exit 1
  }

  if (params.watch == true) {
    epochMap = [:]
    for (i in ["mapping","bamMapping","pairing","aggregate"]) {
      if (file(params."${i}".toString()).exists()){ 
	epochMap[file(params."${i}").toRealPath()] = 0 
      }
    }
    touchInputs(params.chunkSizeLimit, epochMap)
  }

  if (doWF_AggregateFromResult){
    aggregateFromResult(runAggregate, multiqcWesConfig, multiqcWgsConfig, multiqcTempoLogo)
  }
  else{
    //Begin executing modules for the run.
    validate_wf()
    inputMapping = validate_wf.out.inputMapping
    inputPairing = validate_wf.out.inputPairing

    if (doWF_align)
    {
      alignment_wf(inputMapping)
    }

    //Handle input bams as coming originally from bams, or from an alignment this run.
    if (params.bamMapping) {
      inputBam = inputMapping
      if (doWF_QC){
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

    if(doWF_germSV)
    {
      germlineSV_wf(bams)
    }

    if(doWF_facets)
    {
      facets_wf(bamFiles)
    }

    if(doWF_germSNV)
    {
      germlineSNV_wf(bams, bamsTumor, scatter_wf.out.mergedIList, facets_wf.out.facetsForMafAnno)
    }

    if(doWF_loh)
    {
      loh_wf(bams, bamFiles, facets_wf.out.facetsPurity)
    }

    if(doWF_SNV)
    {
      snv_wf(bamFiles, scatter_wf.out.mergedIList, manta_wf.out.mantaToStrelka, loh_wf.out.hlaOutput, facets_wf.out.facetsForMafAnno)
    }

    if(doWF_SV)
    {
      sv_wf(bamFiles, manta_wf.out.manta4Combine)
      if (doWF_SNV && params.assayType == "genome"){
        ascat_wf(bamFiles)
        hrdetect_wf(ascat_wf.out.ascatCNV, snv_wf.out.mafFile, sv_wf.out.dellyMantaCombinedBedpe)
      }
    }

    if(doWF_QC)
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

    if(doWF_QC && params.pairing)
    {
      samplePairingQC_wf(inputBam, inputPairing, runConpairAll)

      FacetsQC4SomaticMultiQC = doWF_facets ? facets_wf.out.FacetsQC4SomaticMultiQC : inputPairing.map{ t_id, n_id -> [t_id, n_id,"",""]}

      qualimap4PairedTN = inputPairing
        .combine(sampleQC_wf.out.qualimap4Process)
        .branch{ idTumor, idNormal, idSample, qualimapDir ->
          tumor: idSample == idTumor
          normal: idSample == idNormal
        }

      qualimap4PairedTN.tumor
        .combine(qualimap4PairedTN.normal, by:[0,1])
        .map{ idTumor,idNormal, idSample1, qualimapDir1, idSample2, qualimapDir2 ->
          [idTumor,idNormal,qualimapDir1,qualimapDir2]
        }.set{qualimap4SomaticMultiQC}


      samplePairingQC_wf.out.conpairOutput
        .map{ placeHolder, idTumor, idNormal, conpairFiles -> [idTumor, idNormal, conpairFiles]}
        .join(qualimap4SomaticMultiQC, by:[0,1])
        .join(FacetsQC4SomaticMultiQC, by:[0,1])
        .set{ somaticMultiQCinput }
        
      somaticMultiQC_wf(somaticMultiQCinput)
    }

    if(doWF_AggregateFromProcess)
    {
      aggregateFromProcess(
        inputPairing,
	runAggregate,
        doWF_facets ? facets_wf : false,
        doWF_SV ? sv_wf : false,
        doWF_SNV ? snv_wf : false,
        doWF_loh ? loh_wf : false,
        doWF_mdParse ? mdParse_wf : false,
        doWF_germSNV ? germlineSNV_wf : false,
        doWF_germSV ? germlineSV_wf : false,
        doWF_QC ? sampleQC_wf : false,
        doWF_QC && params.pairing ? samplePairingQC_wf : false,
        fastPJson,
        multiqcWesConfig, 
        multiqcWgsConfig, 
        multiqcTempoLogo
      )
    }
  }
}
workflow.onComplete {
  file(params.fileTracking).text = ""
  file(outDir).eachFileRecurse{
    file(params.fileTracking).append(it + "\n")
  }
}

