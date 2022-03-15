process GermlineCombineChannel {
  tag {idTumor + "__" + idNormal}

// 3 intermidiate files (plus 3 index files) output for step by step filter check (2 filter steps involved here)
  publishDir "${params.outDir}/germline/${idNormal}/combined_mutations/intermediate_files/", mode: params.publishDirMode, pattern: "*.union.*"
  publishDir "${params.outDir}/germline/${idNormal}/combined_mutations/intermediate_files/", mode: params.publishDirMode, pattern: "*.germline.vcf.gz*"

  input:
    tuple val(idNormal), val(target), val(placeHolder), path(haplotypecallercombinedVcf), path(haplotypecallercombinedVcfIndex), path(strelkaVcf), path(strelkaVcfIndex), val(idTumor), path(bamTumor), path(baiTumor)
    tuple path(genomeFile), path(genomeIndex)
    tuple path(repeatMasker), path(repeatMaskerIndex), path(mapabilityBlacklist), path(mapabilityBlacklistIndex)
    tuple path(gnomadWesVcf), path(gnomadWesVcfIndex), path(gnomadWgsVcf), path(gnomadWgsVcfIndex)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${idTumor}__${idNormal}.germline.vcf"), emit: mutationMergedGermline
    path("${idNormal}.union.vcf.gz")
    path("${idNormal}.union.vcf.gz.tbi")
    path("${idNormal}.union.pass.vcf.gz")
    path("${idNormal}.union.pass.vcf.gz.tbi")
    path("${idTumor}__${idNormal}.germline.vcf.gz")
    path("${idTumor}__${idNormal}.germline.vcf.gz.tbi")

  script:  
  isecDir = "${idNormal}.isec"
  gnomad = gnomadWgsVcf
  if (params.assayType == 'genome') {
    gnomad = gnomadWgsVcf
  }
  else if (params.assayType == 'exome') {
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=HaplotypeCaller,Number=0,Type=Flag,Description=\"Variant was called by HaplotypeCaller\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e '##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="RepeatMasker">' > vcf.rm.header
  echo -e '##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description="EncodeDacMapability">' > vcf.map.header

  bcftools isec \
    --output-type z \
    --prefix ${isecDir} \
    ${haplotypecallercombinedVcf} ${strelkaVcf}

  bcftools annotate \
    --annotations ${isecDir}/0003.vcf.gz \
    --include 'FILTER!=\"PASS\"' \
    --mark-sites \"+Strelka2FILTER\" \
    -k \
    --output-type z \
    --output ${isecDir}/0003.annot.vcf.gz \
    ${isecDir}/0003.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0000.vcf.gz \
    --mark-sites +HaplotypeCaller \
    --output-type z \
    --output ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0002.vcf.gz \
    --mark-sites \"+HaplotypeCaller;Strelka2\" \
    --output-type z \
    --output ${isecDir}/0002.tmp.vcf.gz \
    ${isecDir}/0002.vcf.gz

  tabix --preset vcf ${isecDir}/0002.tmp.vcf.gz
  tabix --preset vcf ${isecDir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${isecDir}/0003.annot.vcf.gz \
    --columns +FORMAT,Strelka2FILTER \
    --output-type z \
    --output ${isecDir}/0002.annot.vcf.gz \
    ${isecDir}/0002.tmp.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0001.vcf.gz \
    --mark-sites +Strelka2 \
    --output-type z \
    --output ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0001.vcf.gz

  tabix --preset vcf ${isecDir}/0000.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0001.annot.vcf.gz
  tabix --preset vcf ${isecDir}/0002.annot.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0001.annot.vcf.gz \
    ${isecDir}/0002.annot.vcf.gz | \
  bcftools sort | \
  bcftools annotate \
    --header-lines vcf.rm.header \
    --annotations ${repeatMasker} \
    --columns CHROM,FROM,TO,RepeatMasker | \
  bcftools annotate \
    --header-lines vcf.map.header \
    --annotations ${mapabilityBlacklist} \
    --columns CHROM,FROM,TO,EncodeDacMapability \
    --output-type z \
    --output ${idNormal}.union.vcf.gz

  tabix --preset vcf ${idNormal}.union.vcf.gz

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type z \
    --output ${idNormal}.union.pass.vcf.gz \
    ${idNormal}.union.vcf.gz

  tabix --preset vcf ${idNormal}.union.pass.vcf.gz

  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    ${idNormal}.union.pass.vcf.gz | \
  bcftools filter \
    --exclude \"${params.germlineVariant.gnomadAf}\" \
    --output-type v \
    --output ${idNormal}.union.gnomad.vcf 

  GetBaseCountsMultiSample \
    --fasta ${genomeFile} \
    --bam ${idTumor}:${bamTumor} \
    --vcf ${idNormal}.union.gnomad.vcf \
    --output ${idTumor}.genotyped.vcf

  bgzip ${idNormal}.union.gnomad.vcf
  bgzip ${idTumor}.genotyped.vcf
  tabix --preset vcf ${idNormal}.union.gnomad.vcf.gz
  tabix --preset vcf ${idTumor}.genotyped.vcf.gz

  bcftools merge \
    --output ${idTumor}__${idNormal}.germline.vcf \
    --output-type v \
    ${idNormal}.union.gnomad.vcf.gz \
    ${idTumor}.genotyped.vcf.gz

  bgzip -c ${idTumor}__${idNormal}.germline.vcf > ${idTumor}__${idNormal}.germline.vcf.gz
  tabix --preset vcf ${idTumor}__${idNormal}.germline.vcf.gz

  """
}