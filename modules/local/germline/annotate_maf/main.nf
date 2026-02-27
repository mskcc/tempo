// Germline Annotate MAF — matches original Tempo GermlineAnnotateMaf
// vcf2maf with MSKCC-CMO center, VEP88, filtering via R script
process GERMLINE_ANNOTATE_MAF {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/vcf2maf:vep88_1.2.7'

    input:
    tuple val(meta), path(vcf_merged)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(fasta_dict)
    path(vep_cache)
    path(isoforms)

    output:
    tuple val(meta), path("${prefix}.germline.maf"), emit: maf_file
    path("${prefix}.germline.unfiltered.maf"), emit: unfiltered_maf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    perl /opt/vcf2maf.pl \\
        --maf-center MSKCC-CMO \\
        --vep-path /usr/bin/vep \\
        --vep-data ${vep_cache} \\
        --vep-forks 4 \\
        --tumor-id ${meta.tumor_id} \\
        --normal-id ${meta.normal_id} \\
        --vcf-tumor-id ${meta.tumor_id} \\
        --vcf-normal-id ${meta.normal_id} \\
        --input-vcf ${vcf_merged} \\
        --ref-fasta ${fasta} \\
        --custom-enst ${isoforms} \\
        --output-maf ${prefix}.germline.raw.maf \\
        --filter-vcf 0

    Rscript --no-init-file /usr/bin/filter-germline-maf.R \\
        --normal-depth ${params.germline_normal_depth ?: 10} \\
        --normal-vaf ${params.germline_normal_vaf ?: 0.2} \\
        --maf-file ${prefix}.germline.raw.maf \\
        --output-prefix ${prefix}.germline

    cp ${prefix}.germline.raw.maf ${prefix}.germline.unfiltered.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: "1.2.7"
    END_VERSIONS
    """

    stub:
    prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.germline.maf ${prefix}.germline.unfiltered.maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: "1.2.7"
    END_VERSIONS
    """
}
