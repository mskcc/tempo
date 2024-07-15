nextflow.enable.dsl=2
include { AggregateIndels }               from '../process/SNV/AggregateIndels' 
include { RunIndelModel }                 from '../process/SNV/RunIndelModel' 

workflow indel_wf {
  take: 
    bamFiles    
    strelkaOut
    platypusOut
    svabaIndelout

  main:
    // pairingQc    = params.pairing

    targetsMap   = params.targetsMap
    referenceMap = params.referenceMap
    // platypusOut.view()
    // strelkaOut.view()
    // svabaIndelout.view()
    combinedChannel = strelkaOut.combine(platypusOut, by: [0,1,2]).combine(svabaIndelout, by: [0,1,2])

    
    // filter to pass
    AggregateIndels(bamFiles,
                    combinedChannel, 
                    Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.repeatMasker,referenceMap.dbsnp]))

    bamTsvChannel = bamFiles.combine(AggregateIndels.out.tsvGroup, by: [0,1,2])
    bamTsvChannel.view()

    RunIndelModel(bamTsvChannel,
                  Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.repeatMasker,referenceMap.py2bitfile])
    )

    emit:
      indelOut  = RunIndelModel.out.indelOut
}

