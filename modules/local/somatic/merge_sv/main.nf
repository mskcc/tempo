process SOMATIC_MERGE_SV {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::delly=0.8.2 bioconda::bcftools=1.9 bioconda::htslib=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/delly-bcftools:0.0.1' :
        'cmopipeline/delly-bcftools:0.0.1' }"

    input:
    tuple val(meta), path(delly_vcf), path(delly_tbi),
          path(manta_vcf), path(manta_tbi)

    output:
    tuple val(meta), path("${prefix}.delly.manta.vcf.gz"), path("${prefix}.delly.manta.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Reorder Manta samples to match Delly (tumor,normal)
    bcftools view \\
        --samples ${meta.tumor_id},${meta.normal_id} \\
        --output-type z \\
        --output-file ${prefix}.manta.swap.vcf.gz \\
        ${manta_vcf}

    tabix --preset vcf ${prefix}.manta.swap.vcf.gz

    # Concatenate Delly and Manta
    bcftools concat \\
        --allow-overlaps \\
        --output-type z \\
        --output ${prefix}.delly.manta.unfiltered.vcf.gz \\
        ${delly_vcf} ${prefix}.manta.swap.vcf.gz

    tabix --preset vcf ${prefix}.delly.manta.unfiltered.vcf.gz

    # Filter PASS variants on canonical chromosomes
    bcftools filter \\
        --include 'FILTER="PASS"' \\
        --regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,MT,X,Y \\
        ${prefix}.delly.manta.unfiltered.vcf.gz | \\
    bcftools sort \\
        --output-type z \\
        --output-file ${prefix}.delly.manta.vcf.gz

    tabix --preset vcf ${prefix}.delly.manta.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.delly.manta.vcf.gz
    touch ${prefix}.delly.manta.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.20
    END_VERSIONS
    """
}
