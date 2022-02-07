include { GermlineDellyCall }                 from '../GermSV/GermlineDellyCall' 
include { GermlineRunManta }                  from '../GermSV/GermlineRunManta' 
include { GermlineMergeDellyAndManta }        from '../GermSV/GermlineMergeDellyAndManta'

workflow germlineSV_wf
{
  take:
    bams

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypesGermline }

    GermlineDellyCall(svTypesGermline,
                      bams,
                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions]))

    GermlineRunManta(bams,
                     Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                     Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex]))

    GermlineDellyCall.out.dellyFilter4CombineGermline.groupTuple(by: [0,1], size: 5)
            .combine(GermlineRunManta.out.mantaOutputGermline, by: [0,1])
            .set{ dellyMantaChannelGermline }

    GermlineMergeDellyAndManta(dellyMantaChannelGermline,
                               Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    GermlineSVVcf2Bedpe(
      GermlineMergeDellyAndManta.out.dellyMantaCombinedOutputGermline,
      referenceMap.repeatMasker,
      referenceMap.mapabilityBlacklist,
      referenceMap.annotSVref, 
      referenceMap.spliceSites, 
      workflow.projectDir + "/containers/svtools" 
    )

  emit:
    dellyMantaCombined4AggregateGermline    = GermlineMergeDellyAndManta.out.dellyMantaCombined4AggregateGermline
    dellyMantaCombinedTbi4AggregateGermline = GermlineMergeDellyAndManta.out.dellyMantaCombinedTbi4AggregateGermline
}