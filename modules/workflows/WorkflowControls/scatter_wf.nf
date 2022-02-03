include { CreateScatteredIntervals }   from '../Scatter/CreateScatteredIntervals'

workflow scatter_wf
{
  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    targets4Intervals = Channel.from(targetsMap.keySet())
        .map{ targetId ->
          [ targetId, targetsMap."${targetId}".targetsBedGz, targetsMap."${targetId}".targetsBedGzTbi ]
        }

    CreateScatteredIntervals(Channel.value([referenceMap.genomeFile, 
                                            referenceMap.genomeIndex, 
                                            referenceMap.genomeDict]), 
                                            targets4Intervals)
  emit:
    mergedIList = CreateScatteredIntervals.out.mergedIList
}