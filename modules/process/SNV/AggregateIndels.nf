process AggregateIndels {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/indelstuff", mode: params.publishDirMode, pattern: "*{tsv,gz,gz.tbi}"

  input:
    tuple val(combineKey),val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple val(combineKey),path(strelkaVcf),path(strelkaTbi),path(platypusVcf),path(svabaVcf)
    tuple path(refGenome), path(dbsnp), path(fai)

    

  output:
    path("*.tsv"), emit: tsvGroup

  script:
  outputPrefix = "${idTumor}__${idNormal}"

  """

    # split vcfs
    /somaticseq/somaticseq/vcfModifier/splitVcf.py -infile ${strelkaVcf} \
        -snv strelka2.snvs.vcf \
       -indel strelka2_somatic_ind.vcf

    # filter to pass only
    bcftools view -O z -f PASS ${platypusVcf} > Platypus.pass.vcf.gz 

    bcftools view -O z -f PASS strelka2_somatic_ind.vcf > Strelka.pass.vcf.gz 
    bcftools view -O z -f PASS ${svabaVcf} > Svaba.pass.vcf.gz 


    /somaticseq/somaticseq/somaticseq_parallel.py \
        --output-directory    ./ \
        --genome-reference   ${refGenome} \
        --threads           6 \
        --dbsnp-vcf         ${dbsnp} \
        paired \
        --tumor-bam-file    ${bamTumor} \
        --normal-bam-file   ${bamNormal} \
        --strelka-indel     Strelka.pass.vcf.gz \
        --platypus-vcf        Platypus.pass.vcf.gz \
        --arbitrary-indels   Svaba.pass.vcf.gz

      """
}