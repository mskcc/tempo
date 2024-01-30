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
    hlaOutput = RunPolysolver.out.hlaOutput.map{ ["placeHolder"] + it }

    bamFiles.combine(facetsPurity, by: [0,1,2])
            .combine(hlaOutput, by: [1,2])
            .set{ mergedChannelLOHHLA }

    RunLOHHLA(mergedChannelLOHHLA, 
            Channel.value([referenceMap.hlaFasta, referenceMap.hlaDat]))

  emit:
    hlaOutput            = hlaOutput
    lohhla4Aggregate     = RunLOHHLA.out.lohhla4Aggregate
}
