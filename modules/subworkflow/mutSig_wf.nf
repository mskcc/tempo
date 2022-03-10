include { RunMutationSignatures }      from '../process/MutSig/RunMutationSignatures'

workflow mutSig_wf
{
  take:
    mafFile

  main:
    if (!(params.cosmic in ['v2', 'v3'])) {
        println "ERROR: Possible values of mutational signature reference --cosmic is 'v2', 'v3'"
        exit 1
    }
    RunMutationSignatures(mafFile)

  emit:
    mutSig4MetaDataParser = RunMutationSignatures.out.mutSig4MetaDataParser
}
