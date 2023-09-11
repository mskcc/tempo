include { SomaticRunStrelka2 }         from '../process/SNV/SomaticRunStrelka2' 
include { RunMutect2 }                 from '../process/SNV/RunMutect2'
include { SomaticCombineMutect2Vcf }   from '../process/SNV/SomaticCombineMutect2Vcf' 
include { SomaticCombineChannel }      from '../process/SNV/SomaticCombineChannel' 
include { SomaticAnnotateMaf }         from '../process/SNV/SomaticAnnotateMaf' 
include { RunNeoantigen }              from '../process/SNV/RunNeoantigen' 
include { SomaticFacetsAnnotation }    from '../process/SNV/SomaticFacetsAnnotation' 
include { RunPlatypus }               from '../process/SNV/RunPlatypus' 

workflow snv_wf
{
  take: 
    bamFiles
    mergedIList
    mantaToStrelka
    hlaOutput
    facetsForMafAnno

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    bamFiles.combine(mergedIList, by: 2)
      .map{
        item ->
          def idTumor = item[1]
          def idNormal = item[2]
          def target = item[0]
          def tumorBam = item[3]
          def normalBam = item[4]
          def tumorBai = item[5]
          def normalBai = item[6]
          def intervalBed = item[7]
          def key = idTumor+"__"+idNormal+"@"+target // adding one unique key

          return [ key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed ]
      }.map{ 
        key, idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed -> 
        tuple ( 
            groupKey(key, intervalBed.size()), // adding numbers so that each sample only wait for it's own children processes
            idTumor, idNormal, target, tumorBam, normalBam, tumorBai, normalBai, intervalBed
        )
    }
    .transpose()
    .set{ mergedChannelSomatic }


    RunPlatypus(mergedChannelSomatic, 
          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    RunMutect2(mergedChannelSomatic, 
          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    RunMutect2.out.forMutect2Combine.groupTuple().set{ forMutect2Combine }

    SomaticCombineMutect2Vcf(forMutect2Combine, 
                             Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    bamFiles.combine(mantaToStrelka, by: [0, 1, 2])
        .map{ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi ->
              [idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi]
        }.set{ input4Strelka }

    SomaticRunStrelka2(input4Strelka,
                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    SomaticCombineMutect2Vcf.out.mutect2CombinedVcfOutput.combine(bamFiles, by: [0,1,2]).combine(SomaticRunStrelka2.out.strelka4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }

    SomaticCombineChannel(mutectStrelkaChannel,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.exomePoN, referenceMap.wgsPoN,referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]))

    SomaticAnnotateMaf(SomaticCombineChannel.out.mutationMergedVcf,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
                                        referenceMap.vepCache, referenceMap.isoforms]))


    hlaOutput.combine(SomaticAnnotateMaf.out.mafFile, by: [1,2]).set{ input4Neoantigen }

    RunNeoantigen(input4Neoantigen, Channel.value([referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS]))

    facetsForMafAnno.combine(RunNeoantigen.out.mafFileForMafAnno, by: [0,1,2]).set{ facetsMafFileSomatic }

    SomaticFacetsAnnotation(facetsMafFileSomatic)

  emit:
    strelkaOut            = SomaticRunStrelka2.out.strelka4Combine
    platypusOut           = RunPlatypus.out.platypusOutput
    mafFile               = SomaticAnnotateMaf.out.mafFile
    maf4MetaDataParser    = SomaticFacetsAnnotation.out.maf4MetaDataParser
    NetMhcStats4Aggregate = RunNeoantigen.out.NetMhcStats4Aggregate
    finalMaf4Aggregate    = SomaticFacetsAnnotation.out.finalMaf4Aggregate
}
