include { RunPolysolver }              from '../process/hlaTyping/RunPolysolver'

workflow hlaTyping_wf
{
  take:
    bams

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    RunPolysolver(bams)
    hlaOutput = RunPolysolver.out.hlaOutput.map{ ["placeHolder"] + it }

  emit:
    hlaOutput            = hlaOutput
}
