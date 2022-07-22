include { GermlineDellyCall }                 from '../process/GermSV/GermlineDellyCall' 
include { DellyCombine
            as GermlineDellyCombine }         from '../process/SV/DellyCombine'
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
    GermlineDellyCombine(
      GermlineDellyCall.out.dellyFilter4CombineGermline
        .groupTuple( by: [0,1], size: 5 )
        .map{ normal_id, target, vcf, tbi ->
          [ "", normal_id, target, vcf, tbi ]
        }
      , "germline"
    )

    GermlineRunManta(
        bams,
        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
        Channel.value([referenceMap.svCallingIncludeRegions, referenceMap.svCallingIncludeRegionsIndex])
    )

    GermlineDellyCombine.out
        .map{ tumor_id, normal_id, target, vcf, tbi -> [normal_id, target, vcf, tbi, "delly" ] }
        .mix(GermlineRunManta.out.mantaOutputGermline.map{ it + ["manta"]})
        .groupTuple( by:[0,1], size:2 )
        .set{allSvCallsCombineChannel}

    GermlineMergeSVs(
        allSvCallsCombineChannel,
			  workflow.projectDir + "/containers/bcftools-vt-mergesvvcf"
		)

    GermlineSVVcf2Bedpe(
        GermlineMergeSVs.out.SVsCombinedOutputGermline
    )
    GermlineAnnotateSVBedpe(
        GermlineSVVcf2Bedpe.out.GermlineCombinedUnfilteredBedpe,
        referenceMap.repeatMasker,
        referenceMap.mapabilityBlacklist,
        referenceMap.spliceSites, 
        workflow.projectDir + "/containers/iannotatesv",
        params.genome
    )
    
    GermlineMergeSVs.out.SVsCombinedOutputGermline
        .map{ idNormal, target, vcfFile, tbiFile ->
          ["placeHolder", "noTumor", idNormal, vcfFile, tbiFile]
        }.set{sv4AggregateGermline}

  emit:
    SVAnnotBedpe         = GermlineAnnotateSVBedpe.out.SVAnnotBedpe
    SVAnnotBedpePass     = GermlineAnnotateSVBedpe.out.SVAnnotBedpePass
    sv4AggregateGermline = sv4AggregateGermline
}
