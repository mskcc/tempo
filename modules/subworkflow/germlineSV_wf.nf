include { GermlineDellyCall }                 from '../process/GermSV/GermlineDellyCall' 
include { GermlineRunManta }                  from '../process/GermSV/GermlineRunManta' 
include { GermlineMergeSVs }                  from '../process/GermSV/GermlineMergeSVs'
include { GermlineSVVcf2Bedpe }               from '../process/GermSV/GermlineSVVcf2Bedpe'

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

    GermlineMergeSVs(dellyMantaChannelGermline,
                               Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]),
			       workflow.projectDir + "/containers/bcftools-vt-mergesvvcf"
			       )

    if (workflow.profile != "test" && workflow.profile != "test_singularity") {
      GermlineSVVcf2Bedpe(
      GermlineMergeSVs.out.SVsCombinedOutputGermline,
      referenceMap.repeatMasker,
      referenceMap.mapabilityBlacklist,
      referenceMap.annotSVref, 
      referenceMap.spliceSites, 
      workflow.projectDir + "/containers/svtools" 
    )
    }

  emit:
    sv4AggregateGermline    = GermlineMergeSVs.out.SVsCombinedOutputGermline
}
