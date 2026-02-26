process VCF2MAF {
    tag "$meta.id"
    label 'process_medium'

    container "cmopipeline/vcf2maf:vep88_1.2.7"

    input:
    tuple val(meta), path(vcf)
    path fasta
    path vep_cache

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_id  = meta.tumor_id  ?: meta.id
    def normal_id = meta.normal_id ?: ''
    """
    vcf2maf.pl \\
        --input-vcf ${vcf} \\
        --output-maf ${prefix}.maf \\
        --tumor-id ${tumor_id} \\
        ${normal_id ? "--normal-id ${normal_id}" : ''} \\
        --ref-fasta ${fasta} \\
        --vep-path /usr/local/bin \\
        --vep-data ${vep_cache} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: \$(vcf2maf.pl --help 2>&1 | grep -i version | head -1 | sed 's/.*version //' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: "stub"
    END_VERSIONS
    """
}
