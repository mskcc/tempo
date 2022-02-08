process DoFacets {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "*.snp_pileup.gz"
  publishDir "${params.outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "${tag}_OUT.txt"
  publishDir "${params.outDir}/somatic/${tag}/facets/${tag}", mode: params.publishDirMode, pattern: "${outputDir}/*.{Rdata,png,out,seg,txt}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    file(facetsVcf)

  output:
    path("${outfile}"), emit: snpPileupOutput
    path("${outputDir}/*"), emit: FacetsOutput
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("*/*_purity.seg"), path("*/*_hisens.seg"), path("*_OUT.txt"), path("*/*.arm_level.txt"), path("*/*.gene_level.txt"), emit: facets4Aggregate
    tuple val(idTumor), val(idNormal), val(target), path("${outputDir}/*purity.out"), emit: facetsPurity
    tuple val(idTumor), val(idNormal), val(target), path("${outputDir}/*purity.Rdata"), path("${outputDir}/*purity.cncf.txt"), path("${outputDir}/*hisens.cncf.txt"), val("${outputDir}"), emit: facetsForMafAnno
    tuple val(idTumor), val(idNormal), val(target), path("${outputDir}/*.{Rdata,png,out,seg,txt}"), path("${idTumor}__${idNormal}.snp_pileup.gz"), val("${outputDir}"), emit: Facets4FacetsPreview
    tuple val("placeHolder"), val(idTumor), val(idNormal), path("*/*.*_level.txt"), emit: FacetsArmGeneOutput
    tuple val(idTumor), val(idNormal), val(target), path("*/*.qc.txt"), emit: FacetsQC4MetaDataParser
    tuple val(idTumor), val(idNormal), path("*_OUT.txt"), emit: FacetsRunSummary

  script:
  tag = outputFacetsSubdirectory = "${idTumor}__${idNormal}"
  outfile = tag + ".snp_pileup.gz"
  outputDir = "facets${params.facets.R_lib}c${params.facets.cval}pc${params.facets.purity_cval}"

  """
  touch .Rprofile

  export SNP_PILEUP=/usr/bin/snp-pileup

  Rscript /usr/bin/facets-suite/snp-pileup-wrapper.R \
    --pseudo-snps 50 \
    --vcf-file ${facetsVcf} \
    --output-prefix ${tag} \
    --normal-bam ${bamNormal} \
    --tumor-bam ${bamTumor}

  mkdir ${outputDir}

  set +e
  i=1
  seed=\$((${params.facets.seed}-1))
  attemptNumber=0

  while [ \$i -eq 1 ]
  do
  attemptNumber=\$(( attemptNumber + 1 ))
  if [ \$attemptNumber -gt 4 ]; then 
    break
  fi
  seed=\$((seed+i))
  
  Rscript /usr/bin/facets-suite/run-facets-wrapper.R \
    --cval ${params.facets.cval} \
    --snp-window-size ${params.facets.snp_nbhd} \
    --normal-depth ${params.facets.ndepth} \
    --min-nhet ${params.facets.min_nhet} \
    --purity-cval ${params.facets.purity_cval}\
    --purity-min-nhet ${params.facets.purity_min_nhet} \
    --genome ${params.facets.genome} \
    --counts-file ${outfile} \
    --sample-id ${tag} \
    --directory ${outputDir} \
    --facets-lib-path /usr/local/lib/R/site-library \
    --seed \$seed \
    --everything \
    --legacy-output T

  i=\$?
  done
  set -e

  python3 /usr/bin/summarize_project.py \
    -p ${tag} \
    -c ${outputDir}/*cncf.txt \
    -o ${outputDir}/*out \
    -s ${outputDir}/*seg
  """
}
