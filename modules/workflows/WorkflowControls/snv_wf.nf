include { RunStrelka2 }                from '../SNV/RunStrelka2' 
include { RunMutect2 }                 from '../SNV/RunMutect2'
include { CombineMutect2Vcf }          from '../SNV/CombineMutect2Vcf' 
include { CombineChannel }             from '../SNV/CombineChannel' 
include { AnnotateMaf }                from '../SNV/AnnotateMaf' 
include { RunNeoantigen }              from '../SNV/RunNeoantigen' 
include { FacetsAnnotation }           from '../SNV/FacetsAnnotation' 

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

    RunMutect2(mergedChannelSomatic, 
          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    RunMutect2.out.forMutect2Combine.groupTuple().set{ forMutect2Combine }

    CombineMutect2Vcf(forMutect2Combine, 
                             Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    bamFiles.combine(mantaToStrelka, by: [0, 1, 2])
        .map{ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi ->
              [idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal, mantaCSI, mantaCSIi, targetsMap."$target".targetsBedGz, targetsMap."$target".targetsBedGzTbi]
        }.set{ input4Strelka }

    RunStrelka2(input4Strelka,
                      Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    CombineMutect2Vcf.out.mutect2CombinedVcfOutput.combine(bamFiles, by: [0,1,2]).combine(RunStrelka2.out.strelka4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }

    CombineChannel(mutectStrelkaChannel,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.exomePoN, referenceMap.wgsPoN,referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]))

    AnnotateMaf(CombineChannel.out.mutationMergedVcf,
                        Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict,
                                        referenceMap.vepCache, referenceMap.isoforms]))


    hlaOutput.combine(AnnotateMaf.out.mafFile, by: [1,2]).set{ input4Neoantigen }

    RunNeoantigen(input4Neoantigen, Channel.value([referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS]))

    facetsForMafAnno.combine(RunNeoantigen.out.mafFileForMafAnno, by: [0,1,2]).set{ facetsMafFileSomatic }

    FacetsAnnotation(facetsMafFileSomatic)

  emit:
    mafFile               = AnnotateMaf.out.mafFile
    maf4MetaDataParser    = FacetsAnnotation.out.maf4MetaDataParser
    NetMhcStats4Aggregate = RunNeoantigen.out.NetMhcStats4Aggregate
    finalMaf4Aggregate    = FacetsAnnotation.out.finalMaf4Aggregate
}