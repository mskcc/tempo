include { RunLOHHLA }                  from '../process/hlaLoH/RunLOHHLA' 

workflow hlaLoH_wf
{
  take:
    hlaOutput
    bamFiles
    facetsPurity

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    bamFiles.combine(facetsPurity, by: [0,1,2])
            .combine(hlaOutput, by: [1,2])
            .set{ mergedChannelLOHHLA }

    RunLOHHLA(mergedChannelLOHHLA, 
            Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat]))

  emit:
    lohhla4Aggregate     = RunLOHHLA.out.lohhla4Aggregate
}
