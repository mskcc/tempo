process DELLY_COMBINE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::delly=0.8.2 bioconda::bcftools=1.9 bioconda::htslib=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/delly-bcftools:0.0.1' :
        'cmopipeline/delly-bcftools:0.0.1' }"

    input:
    tuple val(meta), path(vcfs), path(tbis)

    output:
    tuple val(meta), path("${prefix}.delly.vcf.gz"), path("${prefix}.delly.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools concat \\
        --allow-overlaps \\
        --output-type z \\
        --output ${prefix}.delly.vcf.gz \\
        ${args} \\
        ${vcfs}

    tabix --preset vcf ${prefix}.delly.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.delly.vcf.gz
    touch ${prefix}.delly.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
