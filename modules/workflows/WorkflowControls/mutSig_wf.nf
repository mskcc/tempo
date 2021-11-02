include { RunMutationSignatures }      from '../MutSig/RunMutationSignatures'

workflow mutSig_wf
{
  take:
    mafFile

  main:
      RunMutationSignatures(mafFile)

  emit:
    mutSig4MetaDataParser = RunMutationSignatures.out.mutSig4MetaDataParser
}