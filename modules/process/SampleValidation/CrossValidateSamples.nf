process CrossValidateSamples {
  input:
    val(inputMapping)
    val(inputPairing)
    
  output:
    val(validSamples), emit: validSamples
    val(invalidSamples), emit: invalidSamples
    val(validPairings), emit: validPairings

  exec:
    //Restructure pairs.
    pairsList = []
    curPair = []
    mod = 0
    for(String pairItem : inputPairing)
    {
      curPair.add(pairItem)
      if(mod == 1)
      {
        pairsList.add(curPair)
        mod = 0
        curPair = []
        continue
      }
      mod = mod + 1
    } 

    //Restructure mappings.
    mappingsList = []
    curMapping = []
    mod = 0
    for(String mapItem : inputMapping)
    {
      curMapping.add(mapItem)
      if(mod == 3)
      {
        mappingsList.add(curMapping)
        mod = 0
        curMapping = []
        continue
      }
      mod = mod + 1
    } 

    //Identify samples where both pairing ids have valid mappings.
    validSamplesList = []
    validPairings    = []
    for(i in 0..pairsList.size()-1)
    {
      if (inputMapping.toString().contains(pairsList[i][0]) && inputMapping.toString().contains(pairsList[i][1]))
      {
        validSamplesList.add(pairsList[i][0])
        validSamplesList.add(pairsList[i][1])
        validPairings.add([pairsList[i][0],pairsList[i][1]])
      }
    }
    
    //If nothing is valid throw and error and stop execution.
    if(validSamplesList.size() == 0)
    {
      println "Error: CrossValidateSamples - No valid samples identified between pairing and mapping files."
      sleep(500)
      exit 1
    }

    //Create valid inputMappings.
    validInputMappings = []
    invalidInputMappings = []
    for(i in 0..mappingsList.size()-1)
    {
      if(validSamplesList.toString().contains(mappingsList[i][0]))
      {
        validInputMappings.add(mappingsList[i])
      }
      else{
        invalidInputMappings.add(mappingsList[i])
      }
    }
    validSamples   = validInputMappings
    invalidSamples = invalidInputMappings

}
