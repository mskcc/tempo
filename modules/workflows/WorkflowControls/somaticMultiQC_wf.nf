include { SomaticRunMultiQC }                  from '../QC/SomaticRunMultiQC'

workflow somaticMultiQC_wf
{
  take:
    somaticMultiQCinput

  main:
    SomaticRunMultiQC(somaticMultiQCinput, Channel.value([params.multiqcWesConfig, params.multiqcWgsConfig, params.multiqcTempoLogo]))
}