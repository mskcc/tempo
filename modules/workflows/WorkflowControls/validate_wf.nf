include { CrossValidateSamples }       from '../SampleValidation/CrossValidateSamples'

workflow validate_wf
{
  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    
    TempoUtils.checkAssayType(params.assayType)
    target_id_list = targetsMap.keySet()
    if (params.watch == false) {
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: true) : file(params.bamMapping, checkIfExists: true)      
      inputMapping = params.mapping ? TempoUtils.extractFastq(mappingFile, params.assayType, target_id_list) : TempoUtils.extractBAM(mappingFile, params.assayType, target_id_list)
    }
    else if (params.watch == true) {
      mappingFile = params.mapping ? file(params.mapping, checkIfExists: false) : file(params.bamMapping, checkIfExists: false)
      inputMapping  = params.mapping ? watchMapping(mappingFile, params.assayType, target_id_list) : watchBamMapping(mappingFile, params.assayType, target_id_list)
    }
    else{}
    if(params.pairing){
      if (params.watch == false) {
        pairingFile = file(params.pairing, checkIfExists: true)
        inputPairing  = TempoUtils.extractPairing(pairingFile)
        TempoUtils.crossValidateTargets(inputMapping, inputPairing)

        CrossValidateSamples(inputMapping.collect(), inputPairing.collect())
        CrossValidateSamples.out.validSamples
          .flatten()
          .collate(4)
          .map { idSample, target, file_pe1, file_pe2 ->
            [idSample, target, file(file_pe1), file(file_pe2)]
          }
          .set { inputMapping }
        CrossValidateSamples.out.validPairings.flatten().collate(2).set{ inputPairing }
      }
      else if (params.watch == true) {
        pairingFile = file(params.pairing, checkIfExists: false)
        inputPairing  = watchPairing(pairingFile)
      }
      else{}
    }
    else {
      inputPairing = Channel.empty()
    }

  emit:
    inputMapping = inputMapping
    inputPairing = inputPairing
}
