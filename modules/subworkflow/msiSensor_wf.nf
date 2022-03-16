include { RunMsiSensor }               from '../process/MSI/RunMsiSensor'

workflow msiSensor_wf
{
  take:
    bamFiles

  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap
    RunMsiSensor(bamFiles,
                Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.genomeDict, referenceMap.msiSensorList]))

  emit:
    msi4MetaDataParser = RunMsiSensor.out.msi4MetaDataParser
}
