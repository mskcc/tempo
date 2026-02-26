process DELLY_CALL_GERMLINE {
    tag "${meta.id}@${svType}"
    label 'process_medium'

    conda "bioconda::delly=0.8.2 bioconda::bcftools=1.9 bioconda::htslib=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/delly-bcftools:0.0.1' :
        'cmopipeline/delly-bcftools:0.0.1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    each svType
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(exclude_bed)

    output:
    tuple val(meta), val(svType),
          path("${prefix}_${svType}.delly.vcf.gz"),
          path("${prefix}_${svType}.delly.vcf.gz.tbi"),
          emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_cmd = exclude_bed ? "--exclude ${exclude_bed}" : ''
    """
    delly call \\
        --svtype ${svType} \\
        --genome ${fasta} \\
        ${exclude_cmd} \\
        --outfile ${prefix}_${svType}.bcf \\
        ${bam}

    delly filter \\
        --filter germline \\
        --outfile ${prefix}_${svType}.filter.bcf \\
        ${prefix}_${svType}.bcf

    bcftools view --output-type z ${prefix}_${svType}.filter.bcf > ${prefix}_${svType}.delly.vcf.gz
    tabix --preset vcf ${prefix}_${svType}.delly.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: \$(delly --version 2>&1 | sed -n '1s/Delly version: *v//p')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_${svType}.delly.vcf.gz
    touch ${prefix}_${svType}.delly.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        delly: 1.2.6
        bcftools: 1.20
    END_VERSIONS
    """
}
