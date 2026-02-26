process GERMLINE_COMBINE_HC_VCF {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/bcftools-vt:1.1.1' :
        'cmopipeline/bcftools-vt:1.1.1' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("${prefix}.haplotypecaller.vcf.gz"), path("${prefix}.haplotypecaller.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools concat \\
        --allow-overlaps \\
        ${vcfs} | \\
    bcftools sort | \\
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --check-ref s \\
        --multiallelics -both | \\
    bcftools norm --rm-dup all \\
        --output-type z \\
        --output ${prefix}.haplotypecaller.vcf.gz

    tabix --preset vcf ${prefix}.haplotypecaller.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.haplotypecaller.vcf.gz
    touch ${prefix}.haplotypecaller.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
