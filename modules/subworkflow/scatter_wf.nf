include { CreateScatteredIntervals }   from '../process/Scatter/CreateScatteredIntervals'

workflow scatter_wf
{
  take:
    targets4Intervals

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    CreateScatteredIntervals(Channel.value([referenceMap.genomeFile, 
                                            referenceMap.genomeIndex, 
                                            referenceMap.genomeDict]), 
                                            targets4Intervals)
  emit:
    mergedIList = CreateScatteredIntervals.out.mergedIList
}
