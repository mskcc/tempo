include { SomaticDellyCall }           from '../SV/SomaticDellyCall' 
include { SomaticRunSvABA }            from '../SV/SomaticRunSvABA' 
include { SomaticMergeDellyAndManta }  from '../SV/SomaticMergeDellyAndManta' 

workflow sv_wf
{
  take:
    bamFiles
    manta4Combine

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }
    SomaticDellyCall(svTypes, bamFiles,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions]))

    // Put manta output and delly output into the same channel so they can be processed together in the group key
    // that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
    SomaticDellyCall.out.dellyFilter4Combine
      .groupTuple(by: [0,1,2], size: 5)
      .combine(manta4Combine, by: [0,1,2])
      .set{ dellyMantaCombineChannel }

    if (params.assayType == "genome") {
      SomaticRunSvaba(
        bamFiles,
        Channel.value([
          referenceMap.genomeFile, 
          referenceMap.genomeIndex,
          referenceMap.genomeDict,
          referenceMap.bwaIndex
        ])
      )
      dellyMantaCombineChannel
        .combine(SomaticRunSvaba.out.SvABA4Combine, by: [0,1,2])
        .set{dellyMantaCombineChannel}
    } else {
      dellyMantaCombineChannel
        .map{ it + [""]}
        .set{dellyMantaCombineChannel}
    }

    // --- Process Delly and Manta VCFs 
    // Merge VCFs, Delly and Manta
    SomaticMergeDellyAndManta(dellyMantaCombineChannel)

  emit:
    dellyMantaCombined4Aggregate    = SomaticMergeDellyAndManta.out.dellyMantaCombined4Aggregate
    dellyMantaCombinedTbi4Aggregate = SomaticMergeDellyAndManta.out.dellyMantaCombinedTbi4Aggregate
}