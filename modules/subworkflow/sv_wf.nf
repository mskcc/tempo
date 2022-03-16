include { SomaticRunSvABA }            from '../process/SV/SomaticRunSvABA' 
include { brass_wf }                   from './brass_wf' addParams(referenceMap: params.referenceMap)
include { SomaticDellyCall }           from '../process/SV/SomaticDellyCall' 
include { SomaticMergeSVs }            from '../process/SV/SomaticMergeSVs' 
include { SomaticSVVcf2Bedpe }         from '../process/SV/SomaticSVVcf2Bedpe'

workflow sv_wf
{
  take:
    bamFiles
    manta4Combine
    sampleStatistics

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }
    SomaticDellyCall(svTypes, bamFiles,
                     Channel.value([
                       referenceMap.genomeFile, 
                       referenceMap.genomeIndex, 
                       referenceMap.svCallingExcludeRegions
                       ])
                    )

    // Put manta output and delly output into the same channel so they can be processed together in the group key
    // that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
    SomaticDellyCall.out.dellyFilter4Combine
      .groupTuple(by: [0,1,2], size: 5)
      .combine(manta4Combine, by: [0,1,2])
      .set{ dellyMantaCombineChannel }

    if (params.assayType == "genome") {
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
      
      dellyMantaCombineChannel
        .combine(SomaticRunSvABA.out.SvABA4Combine, by: [0,1,2])
        .combine(brass_wf.out.BRASS4Combine, by: [0,1,2])
        .set{allSvCallsCombineChannel}
    } else {
      dellyMantaCombineChannel
        .map{ it + ["","","",""]}
        .set{allSvCallsCombineChannel}
    }

    // --- Process Delly and Manta VCFs 
    // Merge VCFs, Delly and Manta
    SomaticMergeSVs(
      allSvCallsCombineChannel,
      workflow.projectDir + "/containers/bcftools-vt-mergesvvcf"
    )

    SomaticSVVcf2Bedpe(
      SomaticMergeSVs.out.SVCallsCombinedVcf,
      referenceMap.repeatMasker,
      referenceMap.mapabilityBlacklist,
      referenceMap.annotSVref, 
      referenceMap.spliceSites, 
      workflow.projectDir + "/containers/svtools" 
    )

  emit:
    SVCombinedBedpe         = SomaticSVVcf2Bedpe.out.SVCombinedBedpe
    SVCombinedBedpePass     = SomaticSVVcf2Bedpe.out.SVCombinedBedpePass
    sv4Aggregate            = SomaticMergeSVs.out.SVCallsCombinedVcf
}
