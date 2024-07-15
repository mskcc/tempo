process AggregateIndels {
  tag "${idTumor + "__" + idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/IndelModelWF", mode: params.publishDirMode

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
    tuple val(idTumor), val(idNormal), val(target), path(strelkaVcf), path(strelkaTbi), path(platypusVcf), path(svabaVcf)
    tuple path(genomeFile), path(genomeIndex), path(genomeDict), path(dbsnp)

    

  output:
    tuple val(idTumor), val(idNormal), val(target), path("Ensemble.sINDEL.tsv"), emit: tsvGroup

  script:
  outputPrefix = "${idTumor}__${idNormal}"

  """
    ls
    # split vcfs
    bgzip -d -c ${strelkaVcf} > Strelka2.vcf
    /opt/somaticseq/somaticseq/vcfModifier/splitVcf.py -infile Strelka2.vcf \
      -snv strelka2.snvs.vcf \
      -indel strelka2_somatic_ind.vcf


    # filter to pass only
    echo "Platypus filter"

    
    bcftools view -O z -f PASS ${platypusVcf} > Platypus.pass.vcf.gz 

    tabix -p vcf  Platypus.pass.vcf.gz 
    
    echo "Strelka filter"
    bcftools view -O z -f PASS strelka2_somatic_ind.vcf > Strelka.pass.vcf.gz 
    tabix -p vcf Strelka.pass.vcf.gz 
    
    
    echo "Svaba filter"
    bcftools view -O z -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,x,X" Platypus.pass.vcf.gz  > Platypus.pass.chrs.vcf.gz
    tabix -p vcf  Platypus.pass.chrs.vcf.gz 

    bcftools view -O z -f PASS ${svabaVcf} > Svaba.pass.vcf.gz 

    bgzip -d -c Svaba.pass.vcf.gz  > Svaba.pass.vcf
    sed -i 's/MAPQ/MQ/g' Svaba.pass.vcf

    bgzip -c Svaba.pass.vcf > Svaba.pass.mq.vcf.gz
    tabix -p vcf Svaba.pass.mq.vcf.gz


    /opt/somaticseq/somaticseq/somaticseq_parallel.py \
        --output-directory    ./ \
        --genome-reference   ${genomeFile} \
        --threads           6 \
        --dbsnp-vcf         ${dbsnp} \
        paired \
        --tumor-bam-file    ${bamTumor} \
        --normal-bam-file   ${bamNormal} \
        --strelka-indel     Strelka.pass.vcf.gz \
        --platypus-vcf      Platypus.pass.chrs.vcf.gz  \
        --arbitrary-indels  Svaba.pass.mq.vcf.gz
      """
}