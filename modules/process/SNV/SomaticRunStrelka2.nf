process SomaticRunStrelka2 {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/strelka2", mode: params.publishDirMode, pattern: "*.vcf.{gz,gz.tbi}"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(mantaCSI), path(mantaCSIi), path(targets), path(targetsIndex)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict) 

  output:
    tuple val(idTumor), val(idNormal), val(target), path('*strelka2.vcf.gz'), path('*strelka2.vcf.gz.tbi'), emit: strelka4Combine
    tuple val(combineKey), path('*strelka2.vcf.gz'), path('*strelka2.vcf.gz.tbi'), emit: strelka4IndelCombine
    tuple path('*strelka2.vcf.gz'), path('*strelka2.vcf.gz.tbi'), emit: strelkaOutput

  script:
  options = ""
  intervals = targets
  if (params.assayType == "exome") {
    options = "--exome"
  }
  outputPrefix = "${idTumor}__${idNormal}"
  outfile = "${outputPrefix}.strelka2.vcf.gz"
  combineKey = "${idTumor}__${idNormal}_${target}"

  """
  configureStrelkaSomaticWorkflow.py \
    ${options} \
    --reportEVSFeatures \
    --callRegions ${intervals} \
    --referenceFasta ${genomeFile} \
    --indelCandidates ${mantaCSI} \
    --tumorBam ${bamTumor} \
    --normalBam ${bamNormal} \
    --runDir Strelka

  python Strelka/runWorkflow.py \
    --mode local \
    --jobs ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${outputPrefix}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${outputPrefix}_somatic_snvs.vcf.gz.tbi

  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz Strelka_${outputPrefix}_somatic_snvs.vcf.gz | \
  bcftools reheader \
    --samples samples.txt | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --output-type z \
    --output ${outfile}

  tabix --preset vcf ${outfile}

  if [ \$(vcf-validator ${outfile} 2>&1 | wc -l ) -gt 0 ] ; then exit 1 ; fi
  """
}
