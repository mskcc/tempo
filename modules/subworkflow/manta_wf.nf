include { SomaticRunManta }            from '../process/SV/SomaticRunManta' 

workflow manta_wf
{
  take: bamFiles

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    SomaticRunManta(bamFiles, 
              Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
              Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex]))
  emit:
    manta4Combine    = SomaticRunManta.out.manta4Combine
    mantaOutput      = SomaticRunManta.out.mantaOutput
    mantaToStrelka   = SomaticRunManta.out.mantaToStrelka
}
