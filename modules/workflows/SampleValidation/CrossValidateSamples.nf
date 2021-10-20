process CrossValidateSamples {
  input:
    val(samplesInMapping)
    val(samplesInPairing)

  output:
    val(isValid), emit: isValid

  exec:
    mappingOnly = samplesInMapping - samplesInPairing
    pairingOnly = samplesInPairing - samplesInMapping
    extraSamples = mappingOnly + pairingOnly
    if(extraSamples == [])
    {
      isValid = true
    }
    else
    {
      println "Mapping: ${mappingOnly}"
      println "Pairing: ${pairingOnly}"
      isValid = false
    }
}
