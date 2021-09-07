process MetaDataParser {
  tag {idTumor + "__" + idNormal}
 
  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/meta_data/", mode: params.publishDirMode, pattern: "*.sample_data.txt"

  input:
    tuple val(idNormal), val(target), val(idTumor), path(purityOut), path(mafFile), path(qcOutput), path(msifile), path(mutSig), val(placeHolder), path(polysolverFile), path(codingBed)
    val(runSomatic)

  output:
    path("*.sample_data.txt"), emit: MetaDataOutput
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("*.sample_data.txt"), emit: MetaData4Aggregate

  when: runSomatic

  script:
  codingRegionsBed = codingBed
  """
  create_metadata_file.py \
    --sampleID ${idTumor}__${idNormal} \
    --tumorID ${idTumor} \
    --normalID ${idNormal} \
    --facetsPurity_out ${purityOut} \
    --facetsQC ${qcOutput} \
    --MSIsensor_output ${msifile} \
    --mutational_signatures_output ${mutSig} \
    --polysolver_output ${polysolverFile} \
    --MAF_input ${mafFile} \
    --coding_baits_BED ${codingRegionsBed}
  
  mv ${idTumor}__${idNormal}_metadata.txt ${idTumor}__${idNormal}.sample_data.txt
  """
}