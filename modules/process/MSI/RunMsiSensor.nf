process RunMsiSensor {
  tag "${idTumor + "__" + idNormal}"
  publishDir "${params.outDir}/somatic/${outputPrefix}/MSI", mode: params.publishDirMode, pattern: "*.tsv", enabled: false
  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict), path(msiSensorList) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.msisensor.tsv"), emit: msi4MetaDataParser

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  """
  msisensor msi \
    -d ${msiSensorList} \
    -t ${bamTumor} \
    -n ${bamNormal} \
    -o ${outputPrefix}.msisensor.tsv
  """
}
