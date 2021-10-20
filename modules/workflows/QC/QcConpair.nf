process QcConpair {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outPrefix}/conpair/", mode: params.publishDirMode

  input:
    tuple val(idTumor), val(idNormal), path(pileupTumor), path(pileupNormal)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict)

  output:
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("${outPrefix}.{concordance,contamination}.txt"), emit: conpairOutput
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("${outPrefix}.concordance.txt"), emit: conpairConcord
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("${outPrefix}.contamination.txt"), emit: conpairContami

  script:
  outPrefix = "${idTumor}__${idNormal}"
  conpairPath = "/usr/bin/conpair"

  markersTxt = ""
  if (params.genome == "GRCh37") {
    markersTxt = "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
  }
  else {
    markersTxt = "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt"
  }

  """
  touch .Rprofile # calls to R inside the python scripts make this necessary to avoid loading user .Rprofile
  
  # Make pairing file
  echo "${idNormal}\t${idTumor}" > pairing.txt

   # Verify concordance
  ${conpairPath}/scripts/verify_concordances.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --normal_homozygous_markers_only \
    --outpre=${outPrefix}

  ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \
    --tumor_pileup=${pileupTumor} \
    --normal_pileup=${pileupNormal} \
    --markers=${markersTxt} \
    --pairing=pairing.txt \
    --outpre=${outPrefix}

  mv ${outPrefix}_concordance.txt ${outPrefix}.concordance.txt
  mv ${outPrefix}_contamination.txt ${outPrefix}.contamination.txt
  """
}