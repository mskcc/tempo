process DELLY_CALL_SOMATIC {
    tag "${meta.id}@${svType}"
    label 'process_medium'

    conda "bioconda::delly=0.8.2 bioconda::bcftools=1.9 bioconda::htslib=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/delly-bcftools:0.0.1' :
        'cmopipeline/delly-bcftools:0.0.1' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
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
    # Step 1: Call SVs for this SV type
    delly call \\
        --svtype ${svType} \\
        --genome ${fasta} \\
        ${exclude_cmd} \\
        --outfile ${prefix}_${svType}.bcf \\
        ${tumor_bam} ${normal_bam}

    # Step 2: Create sample classification file
    printf "${meta.tumor_id}\\ttumor\\n${meta.normal_id}\\tcontrol\\n" > samples.tsv

    # Step 3: Somatic filtering
    delly filter \\
        --filter somatic \\
        -a 0.05 \\
        --samples samples.tsv \\
        --outfile ${prefix}_${svType}.filter.bcf \\
        ${prefix}_${svType}.bcf

    # Step 4: Custom read support filtering
    bcftools view \\
        -s ${meta.tumor_id},${meta.normal_id} \\
        ${prefix}_${svType}.filter.bcf | \\
    bcftools filter \\
        --soft-filter tumor_read_supp -m + \\
        -e "FORMAT/DV[0] < 5 | FORMAT/RV[0] < 2" | \\
    bcftools filter \\
        --soft-filter normal_read_supp -m + \\
        -e "FORMAT/DV[1] > 0 | FORMAT/RV[1] > 0" | \\
    bcftools view --output-type z > \\
        ${prefix}_${svType}.delly.vcf.gz

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
