process GERMLINE_ANNOTATE_MAF {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/vcf2maf:vep88_1.2.7'

    input:
    tuple val(meta), path(vcf_merged)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(fasta_dict)
    path(vep_cache)
    path(isoforms)

    output:
    tuple val(meta), path("*.maf"), emit: maf_file
    path("*.unfiltered.maf"), emit: unfiltered_maf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p vcf2maf_output

    perl /opt/vcf2maf.pl \\
        --input-vcf ${vcf_merged} \\
        --output-maf ${prefix}.unfiltered.maf \\
        --ref-fasta ${fasta} \\
        --vep-path /opt/vep \\
        --vep-data ${vep_cache} \\
        --isoforms ${isoforms} \\
        ${args}

    Rscript /opt/filter-germline-maf.R \\
        --input-maf ${prefix}.unfiltered.maf \\
        --output-maf ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: "1.2.7"
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.maf ${prefix}.unfiltered.maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: "1.2.7"
    END_VERSIONS
    """
}
