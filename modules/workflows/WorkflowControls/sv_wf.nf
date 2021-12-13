include { DellyCall }           from '../SV/DellyCall' 
include { MergeDellyAndManta }  from '../SV/MergeDellyAndManta' 

workflow sv_wf
{
  take:
    bamFiles
    manta4Combine

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypes }
    DellyCall(svTypes, bamFiles,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions]))

    // Put manta output and delly output into the same channel so they can be processed together in the group key
    // that they came in with i.e. (`idTumor`, `idNormal`, and `target`)
    DellyCall.out.dellyFilter4Combine.groupTuple(by: [0,1,2], size: 5).combine(manta4Combine, by: [0,1,2]).set{ dellyMantaCombineChannel }

    // --- Process Delly and Manta VCFs 
    // Merge VCFs, Delly and Manta
    MergeDellyAndManta(dellyMantaCombineChannel)

  emit:
    dellyMantaCombined4Aggregate    = MergeDellyAndManta.out.dellyMantaCombined4Aggregate
    dellyMantaCombinedTbi4Aggregate = MergeDellyAndManta.out.dellyMantaCombinedTbi4Aggregate
}