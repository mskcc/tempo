include { RunManta }            from '../SV/RunManta' 

workflow manta_wf
{
  take: bamFiles

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    RunManta(bamFiles, 
            Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
            Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex]))
  emit:
    manta4Combine    = RunManta.out.manta4Combine
    mantaOutput      = RunManta.out.mantaOutput
    mantaToStrelka   = RunManta.out.mantaToStrelka
}