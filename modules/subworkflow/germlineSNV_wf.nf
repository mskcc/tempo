include { GermlineRunHaplotypecaller }        from '../process/GermSNV/GermlineRunHaplotypecaller'
include { GermlineCombineHaplotypecallerVcf } from '../process/GermSNV/GermlineCombineHaplotypecallerVcf' 
include { GermlineRunStrelka2 }               from '../process/GermSNV/GermlineRunStrelka2' 
include { GermlineCombineChannel }            from '../process/GermSNV/GermlineCombineChannel' 
include { GermlineAnnotateMaf }               from '../process/GermSNV/GermlineAnnotateMaf' 
include { GermlineFacetsAnnotation }          from '../process/GermSNV/GermlineFacetsAnnotation' 

workflow germlineSNV_wf
{
  take:
    bams
    bamsTumor
    mergedIList
    facetsForMafAnno

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    bams.combine(mergedIList, by: 1)
      .map{
        item ->
          def idNormal = item[1]
          def target = item[0]
          def normalBam = item[2]
          def normalBai = item[3]
          def intervalBed = item[4]
          def key = idNormal+"@"+target // adding one unique key
          return [ key, idNormal, target, normalBam, normalBai, intervalBed ]
      }.map{
        key, idNormal, target, normalBam, normalBai, intervalBed ->
        tuple (
            groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
            idNormal, target, normalBam, normalBai, intervalBed
        )
      }
      .transpose()
      .set{ mergedChannelGermline }


    GermlineRunHaplotypecaller(mergedChannelGermline,
                              Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    GermlineRunHaplotypecaller.out.haplotypecaller4Combine.groupTuple().set{ haplotypecaller4Combine }

    GermlineCombineHaplotypecallerVcf(haplotypecaller4Combine,
                                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    bams.map{ idNormal, target, bamNormal, baiNormal -> 
              [idNormal, target, bamNormal, baiNormal, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi]
            }.set{ bamsForStrelkaGermline }

    GermlineRunStrelka2(bamsForStrelkaGermline, 
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    // Join HaploTypeCaller and Strelka outputs, bcftools.
    GermlineCombineHaplotypecallerVcf.out.haplotypecallerCombinedVcfOutput.combine(GermlineRunStrelka2.out.strelkaOutputGermline, by: [0,1,2])
            .combine(bamsTumor, by: [1,2])
            .set{ mergedChannelVcfCombine }

    GermlineCombineChannel(mergedChannelVcfCombine,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex,]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex, referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]))

    GermlineAnnotateMaf(GermlineCombineChannel.out.mutationMergedGermline,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
                                       referenceMap.vepCache, referenceMap.isoforms]))

    facetsForMafAnno.combine(GermlineAnnotateMaf.out.mafFileGermline, by: [0,1,2])
      .set{ facetsMafFileGermline }

    GermlineFacetsAnnotation(facetsMafFileGermline)
    
  emit:
    snv4AggregateGermline = GermlineFacetsAnnotation.out.mafFile4AggregateGermline      
}
