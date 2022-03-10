include { RunPolysolver }              from '../process/LoH/RunPolysolver'
include { RunLOHHLA }                  from '../process/LoH/RunLOHHLA' 

workflow loh_wf
{
  take:
    bams
    bamFiles
    facetsPurity

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    RunPolysolver(bams)

    bamFiles.combine(facetsPurity, by: [0,1,2])
            .combine(RunPolysolver.out.hlaOutput, by: [1,2])
            .set{ mergedChannelLOHHLA }

    RunLOHHLA(mergedChannelLOHHLA, 
            Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat]))

  emit:
    hlaOutput            = RunPolysolver.out.hlaOutput
    lohhla4Aggregate     = RunLOHHLA.out.lohhla4Aggregate
}
