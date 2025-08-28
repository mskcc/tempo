process MetaDataParser {
  tag "${idTumor + "__" + idNormal}"
 
  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/meta_data/", mode: params.publishDirMode, pattern: "*.sample_data.txt"

  input:
    tuple val(idNormal), val(target), val(idTumor), path(purityOut), path(mafFile), path(qcOutput), path(msifile), path(mutSig), val(placeHolder), path(polysolverFile), path(codingBed)

  output:
    path("*.sample_data.txt"), emit: MetaDataOutput
    tuple val(placeHolder), val(idTumor), val(idNormal), path("*.sample_data.txt"), emit: MetaData4Aggregate

  script:

  facetsPurity_out = purityOut.size() > 0 ?  "--facetsPurity_out ${purityOut}" : ""
  facetsQC = qcOutput.size() > 0 ? "--facetsQC ${qcOutput}" : ""
  MSIsensor_output = msifile.size() > 0 ? "--MSIsensor_output ${msifile}" : ""
  mutational_signatures_output = mutSig.size() > 0 ? "--mutational_signatures_output ${mutSig}" : ""
  polysolver_output  = polysolverFile.size() > 0 ? "--polysolver_output ${polysolverFile}" : ""
  MAF_input = mafFile.size() > 0 ? "--MAF_input ${mafFile}" : ""
  coding_baits_BED = codingBed.size() > 0 ? "--coding_baits_BED ${codingBed}" : ""

  """
  create_metadata_file.py \
    ${facetsPurity_out} \
    ${facetsQC} \
    ${MSIsensor_output} \
    ${mutational_signatures_output} \
    ${polysolver_output} \
    ${MAF_input} \
    ${coding_baits_BED} \
    --sampleID ${idTumor}__${idNormal} \
    --tumorID ${idTumor} \
    --normalID ${idNormal}
  
  mv ${idTumor}__${idNormal}_metadata.txt ${idTumor}__${idNormal}.sample_data.txt
  """
}
