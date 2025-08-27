include { RunNeoantigen }              from '../process/SNV/RunNeoantigen' 

workflow neoantigen_wf
{
  take: 
    mafFiles
    hlaOutput

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    hlaOutput.combine(mafFiles, by: [1,2]).set{ input4Neoantigen }

    RunNeoantigen(input4Neoantigen, Channel.value([referenceMap.neoantigenCDNA, referenceMap.neoantigenCDS]))

  emit:
    NetMhcStats4Aggregate = RunNeoantigen.out.NetMhcStats4Aggregate
}
