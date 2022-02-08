include { GermlineDellyCall }                 from '../process/GermSV/GermlineDellyCall' 
include { GermlineRunManta }                  from '../process/GermSV/GermlineRunManta' 
include { GermlineMergeDellyAndManta }        from '../process/GermSV/GermlineMergeDellyAndManta'

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

  emit:
    sv4AggregateGermline    = GermlineMergeDellyAndManta.out.sv4AggregateGermline
}
