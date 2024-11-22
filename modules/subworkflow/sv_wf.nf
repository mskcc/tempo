include { SomaticDellyCall }           from '../process/SV/SomaticDellyCall' 
include { DellyCombine
            as SomaticDellyCombine }   from '../process/SV/DellyCombine'
include { SomaticRunSvABA }            from '../process/SV/SomaticRunSvABA' 
include { brass_wf }                   from './brass_wf' addParams(referenceMap: params.referenceMap)
include { SomaticMergeSVs }            from '../process/SV/SomaticMergeSVs' 
include { SomaticSVVcf2Bedpe }         from '../process/SV/SomaticSVVcf2Bedpe'
include { SomaticAnnotateSVBedpe }     from '../process/SV/SomaticAnnotateSVBedpe'
include { SomaticRunClusterSV }        from '../process/SV/SomaticRunClusterSV'
include { RunSVSignatures }            from '../process/HRDetect/RunSVSignatures'
include { SomaticRunSVCircos }         from '../process/SV/SomaticRunSVCircos'

workflow sv_wf
{
  take:
    bamFiles
    manta4Combine
    sampleStatistics
    cnvCalls

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

    SomaticRunSvABA(
      bamFiles
        .map{ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal -> 
	  [ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal ] + [targetsMap."$target".targetsBed] 
	},
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.bwaIndex
    )

    if (params.assayType == "genome" && workflow.profile != "test") {
      brass_wf(
        bamFiles, 
        sampleStatistics // from facets (default) or ascat
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
	.mix(SomaticRunSvABA.out.SvABA4Combine.map{ it + ["svaba"]})
        .groupTuple( by:[0,1,2], size:3 )
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

    if (params.assayType == "genome") {
      RunSVSignatures(
        SomaticAnnotateSVBedpe.out.SVAnnotBedpePass,
        workflow.projectDir + "/containers/signaturetoolslib/sv_signatures_wrapper.R"
      )
      SVSignatures = RunSVSignatures.out
      SomaticRunClusterSV( SomaticAnnotateSVBedpe.out.SVAnnotBedpePass )
      sv4Aggregate = SomaticRunClusterSV.out.Bedpe4Aggregate.map{ ["placeHolder"] + it }
    } else {
      SVSignatures = Channel.empty()
      sv4Aggregate = SomaticAnnotateSVBedpe.out.SVAnnotBedpe4Aggregate.map{ ["placeHolder"] + it }
    }

    if (cnvCalls){
      SomaticRunSVCircos(
        SomaticAnnotateSVBedpe.out.SVAnnotBedpePass
        .combine(cnvCalls, by: [0,1,2]),
        workflow.projectDir + "/containers/biocircos/generate_biocircos.R", 
	workflow.projectDir + "/containers/biocircos/biocircos.Rmd"
      )
    }

  emit:
    SVSignatures         = SVSignatures
    SVAnnotBedpe         = SomaticAnnotateSVBedpe.out.SVAnnotBedpe
    SVAnnotBedpePass     = SomaticAnnotateSVBedpe.out.SVAnnotBedpePass
    sv4Aggregate         = sv4Aggregate
}
