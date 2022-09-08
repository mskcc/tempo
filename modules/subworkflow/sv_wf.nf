include { SomaticDellyCall }           from '../process/SV/SomaticDellyCall' 
include { DellyCombine
            as SomaticDellyCombine }   from '../process/SV/DellyCombine'
include { SomaticRunSvABA }            from '../process/SV/SomaticRunSvABA' 
include { brass_wf }                   from './brass_wf' addParams(referenceMap: params.referenceMap)
include { SomaticMergeSVs }            from '../process/SV/SomaticMergeSVs' 
include { SomaticSVVcf2Bedpe }         from '../process/SV/SomaticSVVcf2Bedpe'
include { SomaticAnnotateSVBedpe }     from '../process/SV/SomaticAnnotateSVBedpe'
include { SomaticRunSVclone }          from '../process/SV/SomaticRunSVclone'
include { SomaticRunClusterSV }        from '../process/SV/SomaticRunClusterSV'
include { RunSVSignatures }            from '../process/HRDetect/RunSVSignatures'

workflow sv_wf
{
  take:
    bamFiles
    manta4Combine
    sampleStatistics
    finalMaf

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }
    SomaticDellyCall(
    	svTypes, 
    	bamFiles,
	Channel.value([
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.svCallingExcludeRegions
	])
    )

    // Put manta output and delly output into the same channel so they can be processed together in the group key
    // that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
  SomaticDellyCombine(
    SomaticDellyCall.out.dellyFilter4Combine
      .groupTuple( by: [0,1,2], size: 5 )
    , "somatic"
  )

    if (params.assayType == "genome" && workflow.profile != "test") {
      SomaticRunSvABA(
        bamFiles,
        referenceMap.genomeFile, 
        referenceMap.genomeIndex,
        referenceMap.genomeDict,
        referenceMap.bwaIndex
      )

      brass_wf(
        bamFiles, 
        sampleStatistics // from ascat
      )
      
      SomaticDellyCombine.out.map{ it + ["delly"]}
        .mix(manta4Combine.map{ it + ["manta"]})
        .mix(SomaticRunSvABA.out.SvABA4Combine.map{ it + ["svaba"]})
        .mix(brass_wf.out.BRASS4Combine.map{ it + ["brass"]})
        .groupTuple( by:[0,1,2], size:4 )
        .set{allSvCallsCombineChannel}
    } else {
      SomaticDellyCombine.out.map{ it + ["delly"]}
        .mix(manta4Combine.map{ it + ["manta"]})
        .groupTuple( by:[0,1,2], size:2 )
        .set{allSvCallsCombineChannel}
    }
    
    // --- Process SV VCFs 
    // Merge VCFs
    SomaticMergeSVs(
      allSvCallsCombineChannel,
      workflow.projectDir + "/containers/bcftools-vt-mergesvvcf"
    )

    // Convert VCF to Bedpe
    SomaticSVVcf2Bedpe(
      SomaticMergeSVs.out.SVCallsCombinedVcf
    )

    // Annotate Bedpe
    SomaticAnnotateSVBedpe(
      SomaticSVVcf2Bedpe.out.SomaticCombinedUnfilteredBedpe,
      referenceMap.repeatMasker,
      referenceMap.mapabilityBlacklist,
      referenceMap.svBlacklistBed,
      referenceMap.svBlacklistBedpe,
      referenceMap.svBlacklistFoldbackBedpe,
      referenceMap.svBlacklistTEBedpe,
      referenceMap.spliceSites,
      workflow.projectDir + "/containers/iannotatesv",
      params.genome
    )

    RunSVSignatures(
      SomaticAnnotateSVBedpe.out.SVAnnotBedpePass,
      workflow.projectDir + "/containers/signaturetoolslib/sv_signatures_wrapper.R"
    )

    SomaticRunSVclone(
      bamFiles
        .combine(SomaticAnnotateSVBedpe.out.SVAnnotBedpePass,by: [0,1,2])
        .combine(finalMaf, by:[0,1,2])
        .combine(sampleStatistics, by: [0,1,2]),
      workflow.projectDir + "/containers/svclone/prepare_svclone_inputs.py"
    )

    SomaticRunClusterSV( SomaticAnnotateSVBedpe.out.SVAnnotBedpePass )

  emit:
    SVAnnotBedpe         = SomaticAnnotateSVBedpe.out.SVAnnotBedpe
    SVAnnotBedpePass     = SomaticAnnotateSVBedpe.out.SVAnnotBedpePass
    sv4Aggregate         = SomaticMergeSVs.out.SVCallsCombinedVcf
}
