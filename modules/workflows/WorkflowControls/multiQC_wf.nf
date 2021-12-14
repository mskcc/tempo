include { RunMultiQC }                  from '../QC/RunMultiQC'

workflow multiQC_wf
{
  take:
    multiQCinput

  main:
    RunMultiQC(multiQCinput, Channel.value([params.multiqcWesConfig, params.multiqcWgsConfig, params.multiqcTempoLogo]))
}