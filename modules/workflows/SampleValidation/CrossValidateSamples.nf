process CrossValidateSamples {
  input:
    val(samplesInMapping)
    val(samplesInPairing)
    val(inputMapping)
    val(inputPairing)

  output:
    val(isValid), emit: isValid
    val(mappingList), emit: validSamples
    val(pairingList), emit: pairingSamples

  exec:
    mappingOnly = samplesInMapping - samplesInPairing
    pairingOnly = samplesInPairing - samplesInMapping
    extraSamples = mappingOnly + pairingOnly
    setDiff = (samplesInMapping + samplesInPairing ) - samplesInMapping.intersect(samplesInPairing)

    //println inputMapping[0]

    if(setDiff.size == 0)
    {
      mappingList = inputMapping
      pairingList = inputPairing
    }
    else{
      validList  = []
      validPairs = []
      for(item in inputMapping)
      {
        if(!setDiff.contains(item[0]))
        {
          validList.add(item)
        }
      }
      for(item in inputPairing)
      {
        if(!setDiff.contains(item[0]) || !setDiff.contains(item[1]))
        {
          validPairs.add(item)
        }
      }
      mappingList = validList
      pairingList = validPairs
    }

    if(extraSamples == [])
    {
      isValid = true
    }
    else
    {
      //println "Mapping: ${mappingOnly}"
      //println "Pairing: ${pairingOnly}"
      isValid = false
    }
}
