include { QcPileup }                           from '../QC/QcPileup'
include { QcConpair }                          from '../QC/QcConpair'

workflow samplePairingQC_wf
{
  take:
    inputChannel
    inputPairing
    runConpairAll

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    QcPileup(inputChannel, Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    QcPileup.out.pileupOutput.combine(inputPairing)
            .filter { item ->
              def idSample = item[0]
              def samplePileup = item[1]
              def idTumor = item[2]
              def idNormal = item[3]
              idSample == idTumor
            }.map { item ->
              def idTumor = item[2]
              def idNormal = item[3]
              def tumorPileup = item[1]
              return [ idTumor, idNormal, tumorPileup ]
            }
            .unique()
            .set{ pileupT }

    QcPileup.out.pileupOutput.combine(inputPairing)
            .filter { item ->
              def idSample = item[0]
              def samplePileup = item[1]
              def idTumor = item[2]
              def idNormal = item[3]
              idSample == idNormal
            }.map { item ->
              def idTumor = item[2]
              def idNormal = item[3]
              def normalPileup = item[1]
              return [ idTumor, idNormal, normalPileup ]
            }
            .unique()
            .set{ pileupN }

    pileupT.combine(pileupN, by: [0, 1]).unique().set{ pileupConpair }

    QcConpair(pileupConpair, Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))

    if(runConpairAll){
      pileupT.combine(pileupN).unique().set{ pileupConpairAll }

      QcConpairAll(pileupConpairAll,
                    Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict]))
    }

  emit:
    // -- Run based on QcConpairAll channels or the single QcConpair channels
    conpair4Aggregate = (!runConpairAll ? QcConpair.out.conpair4Aggregate : QcConpairAll.out.conpairAll4Aggregate)
    conpairOutput = QcConpair.out.conpairOutput
}
