// Germline hard filtering — matches original Tempo GermlineRunHaplotypecaller
// Separates SNPs and INDELs, applies GATK best-practice VariantFiltration, then merges
process GERMLINE_HARD_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::gatk4=4.4.0.0 bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'docker.io/broadinstitute/gatk:4.4.0.0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("${prefix}.haplotypecaller.filtered.vcf.gz"), path("${prefix}.haplotypecaller.filtered.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gatk SelectVariants \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --select-type-to-include SNP \\
        --output ${prefix}.snps.vcf.gz

    gatk SelectVariants \\
        --reference ${fasta} \\
        --variant ${vcf} \\
        --select-type-to-include INDEL \\
        --output ${prefix}.indels.vcf.gz

    gatk VariantFiltration \\
        --reference ${fasta} \\
        --variant ${prefix}.snps.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
        --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
        --filter-expression "FS > 60.0" --filter-name "FS60" \\
        --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        --output ${prefix}.snps.filtered.vcf.gz

    gatk VariantFiltration \\
        --reference ${fasta} \\
        --variant ${prefix}.indels.vcf.gz \\
        --filter-expression "QD < 2.0" --filter-name "QD2" \\
        --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \\
        --filter-expression "FS > 200.0" --filter-name "FS200" \\
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        --output ${prefix}.indels.filtered.vcf.gz

    bcftools concat \\
        --allow-overlaps \\
        ${prefix}.snps.filtered.vcf.gz ${prefix}.indels.filtered.vcf.gz | \\
    bcftools sort \\
        --output-type z \\
        --output ${prefix}.haplotypecaller.filtered.vcf.gz

    tabix --preset vcf ${prefix}.haplotypecaller.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.haplotypecaller.filtered.vcf.gz
    touch ${prefix}.haplotypecaller.filtered.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.4.0.0
        bcftools: 1.20
    END_VERSIONS
    """
}
