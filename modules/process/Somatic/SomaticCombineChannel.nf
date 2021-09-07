process SomaticCombineChannel {
  tag {idTumor + "__" + idNormal}

  // 3 intermidiate files (plus 3 index files) output for step by step filter check (2 filter steps involved here)
  publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/combined_mutations/intermediate_files/", mode: params.publishDirMode, pattern: "*.union.annot.*"

  input:
    tuple val(idTumor), val(idNormal), val(target), path(mutectCombinedVcf), path(mutectCombinedVcfIndex), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal), path(strelkaVcf), path(strelkaVcfIndex)
    tuple path(genomeFile), path(genomeIndex)
    tuple path(repeatMasker), path(repeatMaskerIndex), path(mapabilityBlacklist), path(mapabilityBlacklistIndex)
    tuple path(exomePoN), path(wgsPoN), path(exomePoNIndex), path(wgsPoNIndex)
    tuple path(gnomadWesVcf), path(gnomadWesVcfIndex), path(gnomadWgsVcf), path(gnomadWgsVcfIndex) 
    val(tools)
    val(runSomatic)
    
  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.pass.vcf"), emit: mutationMergedVcf
    path("${idTumor}__${idNormal}.union.annot.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.vcf.gz.tbi")
    path("${idTumor}__${idNormal}.union.annot.filter.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.filter.vcf.gz.tbi")
    path("${idTumor}__${idNormal}.union.annot.filter.pass.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.filter.pass.vcf.gz.tbi")

  when: tools.containsAll(["manta", "strelka2", "mutect2"]) && runSomatic
  
  script:
  outputPrefix = "${idTumor}__${idNormal}"
  isecDir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWgsVcf
  if (target == "wgs") {
    pon = wgsPoN
    gnomad = gnomadWgsVcf
  }
  else {
    pon = exomePoN
    gnomad = gnomadWesVcf
  }
  """
  echo -e "##INFO=<ID=MuTect2,Number=0,Type=Flag,Description=\"Variant was called by MuTect2\">" > vcf.header
  echo -e "##INFO=<ID=Strelka2,Number=0,Type=Flag,Description=\"Variant was called by Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=Strelka2FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in Strelka2\">" >> vcf.header
  echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > vcf.rm.header
  echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > vcf.map.header
  echo -e "##INFO=<ID=PoN,Number=1,Type=Integer,Description=\"Count in panel of normals\">" > vcf.pon.header
  echo -e "##FORMAT=<ID=alt_count_raw,Number=1,Type=Integer,Description=\"Raw alternate allele depth\">" > vcf.ad_n.header
  echo -e "##FORMAT=<ID=alt_count_raw_fwd,Number=1,Type=Integer,Description=\"Raw alternate allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=alt_count_raw_rev,Number=1,Type=Integer,Description=\"Raw alternate allele depth on reverse strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw,Number=1,Type=Integer,Description=\"Raw reference allele depth\">" > vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw_fwd,Number=1,Type=Integer,Description=\"Raw reference allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=ref_count_raw_rev,Number=1,Type=Integer,Description=\"Raw reference allele depth on reverse strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw,Number=1,Type=Integer,Description=\"Raw total allele depth\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw_fwd,Number=1,Type=Integer,Description=\"Raw total allele depth on forward strand\">" >> vcf.ad_n.header
  echo -e "##FORMAT=<ID=depth_raw_rev,Number=1,Type=Integer,Description=\"Raw total allele depth on reverse strand\">" >> vcf.ad_n.header

  # Get set differences of variant calls:
  # 0000: MuTect2 only
  # 0001: Strelka2 only
  # 0002: MuTect2 calls shared by Strelka2
  # 0003: Strelka2 calls shared by MuTect2

  bcftools isec \
    --output-type z \
    --prefix ${isecDir} \
    ${mutectCombinedVcf} ${strelkaVcf}

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
    --mark-sites +MuTect2 \
    --output-type z \
    --output ${isecDir}/0000.annot.vcf.gz \
    ${isecDir}/0000.vcf.gz

  bcftools annotate \
    --header-lines vcf.header \
    --annotations ${isecDir}/0002.vcf.gz \
    --mark-sites \"+MuTect2;Strelka2\" \
    --output-type z \
    --output ${isecDir}/0002.tmp.vcf.gz \
    ${isecDir}/0002.vcf.gz

  tabix --preset vcf ${isecDir}/0002.tmp.vcf.gz
  tabix --preset vcf ${isecDir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${isecDir}/0003.annot.vcf.gz \
    --columns +INFO,+FORMAT,Strelka2FILTER \
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

  # Concatenate the different sets, annotate with blacklists
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
    --output ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.vcf.gz

  # Add gnomAD annotation
  bcftools annotate \
    --annotations ${gnomad} \
    --columns INFO \
    --output-type z \
    --output ${idTumor}.union.gnomad.vcf.gz \
    ${idTumor}.union.vcf.gz

  tabix --preset vcf ${idTumor}.union.gnomad.vcf.gz

  # Add PoN annotation and flanking sequence
  bcftools annotate \
    --header-lines vcf.pon.header \
    --annotations ${pon} \
    --columns PoN:=AC_Het \
    ${idTumor}.union.gnomad.vcf.gz | \
    vt annotate_indels \
    -r ${genomeFile} \
    -o ${idTumor}.union.annot.vcf -

  mv ${idTumor}.union.annot.vcf ${outputPrefix}.union.annot.vcf

  # Do custom filter annotation, then filter variants
  filter-vcf.py ${outputPrefix}.union.annot.vcf

  bcftools filter \
    --include 'FILTER=\"PASS\"' \
    --output-type v \
    --output ${outputPrefix}.vcf \
    ${outputPrefix}.union.annot.filter.vcf

  # Add normal read count, using all reads
  GetBaseCountsMultiSample \
    --thread ${task.cpus} \
    --maq 0 \
    --fasta ${genomeFile} \
    --bam ${idTumor}:${bamTumor} \
    --bam ${idNormal}:${bamNormal} \
    --vcf ${outputPrefix}.vcf \
    --output ${outputPrefix}.genotyped.vcf 
  
  bgzip ${outputPrefix}.vcf
  bgzip ${outputPrefix}.genotyped.vcf
  tabix --preset vcf ${outputPrefix}.vcf.gz
  tabix --preset vcf ${outputPrefix}.genotyped.vcf.gz

  bcftools annotate \
    --annotations ${outputPrefix}.genotyped.vcf.gz \
    --header-lines vcf.ad_n.header \
    --columns FORMAT/alt_count_raw:=FORMAT/AD,FORMAT/ref_count_raw:=FORMAT/RD,FORMAT/alt_count_raw_fwd:=FORMAT/ADP,FORMAT/ref_count_raw_fwd:=FORMAT/RDP,FORMAT/alt_count_raw_rev:=FORMAT/ADN,FORMAT/ref_count_raw_rev:=FORMAT/RDN,FORMAT/depth_raw:=FORMAT/DP,FORMAT/depth_raw_fwd:=FORMAT/DPP,FORMAT/depth_raw_rev:=FORMAT/DPN \
    --output-type v \
    --output ${outputPrefix}.union.annot.filter.pass.vcf \
    ${outputPrefix}.vcf.gz

  cp ${outputPrefix}.union.annot.filter.pass.vcf ${outputPrefix}.pass.vcf
  bgzip ${outputPrefix}.union.annot.vcf
  bgzip ${outputPrefix}.union.annot.filter.vcf
  bgzip ${outputPrefix}.union.annot.filter.pass.vcf
  tabix --preset vcf ${outputPrefix}.union.annot.vcf.gz
  tabix --preset vcf ${outputPrefix}.union.annot.filter.vcf.gz
  tabix --preset vcf ${outputPrefix}.union.annot.filter.pass.vcf.gz

  """
}