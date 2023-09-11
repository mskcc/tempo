

workflow indel_wf
{
  take: 
    bamFiles    
    strelkaOut
    platypusOut
    svabaIndelout

  main:
    
    combinedChannel = strelkaOut.combine(platypusOut, by:[0,1,2]).combine(svabaIndelout, by:[0,1,2])
    combinedChannel.view()


}
