include { GermlineDellyCall }                 from '../process/GermSV/GermlineDellyCall' 
include { GermlineRunManta }                  from '../process/GermSV/GermlineRunManta' 
include { GermlineMergeSVs }                  from '../process/GermSV/GermlineMergeSVs'
include { GermlineSVVcf2Bedpe }               from '../process/GermSV/GermlineSVVcf2Bedpe'
include { GermlineAnnotateSVBedpe }           from '../process/GermSV/GermlineAnnotateSVBedpe'

workflow germlineSV_wf
{
  take:
    bams

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    Channel.from("DUP", "BND", "DEL", "INS", "INV").set{ svTypesGermline }

    GermlineDellyCall(
        svTypesGermline,
        bams,
        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.svCallingExcludeRegions])
    )

    GermlineRunManta(
        bams,
        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
        Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex])
    )

    GermlineDellyCall.out.dellyFilter4CombineGermline.groupTuple(by: [0,1], size: 5)
        .combine(GermlineRunManta.out.mantaOutputGermline, by: [0,1])
        .set{ dellyMantaChannelGermline }

    GermlineMergeSVs(
        dellyMantaChannelGermline,
			  workflow.projectDir + "/containers/bcftools-vt-mergesvvcf"
		)

    GermlineSVVcf2Bedpe(
        GermlineMergeSVs.out.SVsCombinedOutputGermline,
    )
    GermlineAnnotateSVBedpe(
        GermlineSVVcf2Bedpe.out.GermlineCombinedUnfilteredBedpe,
        referenceMap.repeatMasker,
        referenceMap.mapabilityBlacklist,
        referenceMap.spliceSites, 
        workflow.projectDir + "/containers/iannotatesv",
        params.genome
    )
    

  emit:
    SVAnnotBedpe         = GermlineAnnotateSVBedpe.out.SVAnnotBedpe
    SVAnnotBedpePass     = GermlineAnnotateSVBedpe.out.SVAnnotBedpePass
    sv4AggregateGermline = GermlineMergeSVs.out.SVsCombinedOutputGermline
}
