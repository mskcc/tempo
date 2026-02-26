process STRELKA2_COMBINE_SOMATIC {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::strelka=2.9.10 bioconda::manta=1.5.0 bioconda::bcftools=1.9 bioconda::vt=0.57721"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/strelka2-manta-bcftools-vt:2.0.1' :
        'cmopipeline/strelka2-manta-bcftools-vt:2.0.1' }"

    input:
    tuple val(meta), path(snv_vcf), path(snv_tbi), path(indel_vcf), path(indel_tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.strelka2.vcf.gz"), path("${prefix}.strelka2.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -e 'TUMOR ${meta.tumor_id}\\nNORMAL ${meta.normal_id}' > samples.txt

    bcftools concat \\
        --allow-overlaps \\
        ${indel_vcf} ${snv_vcf} | \\
    bcftools reheader \\
        --samples samples.txt | \\
    bcftools sort | \\
    bcftools norm \\
        --fasta-ref ${fasta} \\
        --check-ref s \\
        --output-type z \\
        --output ${prefix}.strelka2.vcf.gz

    tabix --preset vcf ${prefix}.strelka2.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.strelka2.vcf.gz
    touch ${prefix}.strelka2.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
