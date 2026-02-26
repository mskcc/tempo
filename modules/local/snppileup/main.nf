process SNPPILEUP {
    tag "$meta.id"
    label 'process_medium'

    container "cmopipeline/facets-suite-preview-htstools:0.0.1"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path facets_vcf

    output:
    tuple val(meta), path("*.snp_pileup.gz"), emit: pileup
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    snp-pileup \\
        ${args} \\
        --gzip \\
        ${facets_vcf} \\
        ${prefix}.snp_pileup.gz \\
        ${normal_bam} \\
        ${tumor_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-pileup: \$(snp-pileup --help 2>&1 | head -1 | sed 's/.*version //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.snp_pileup.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-pileup: "stub"
    END_VERSIONS
    """
}
