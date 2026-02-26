process GERMLINE_FACETS_ANNOTATION {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/facets-suite-preview-htstools:0.0.1'

    input:
    tuple val(meta), path(hisens_rdata), path(maf)

    output:
    tuple val(meta), path("*.final.maf"), emit: final_maf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p annotation_output

    # Check if MAF file is empty before processing
    MAF_LINES=\$(wc -l < ${maf})
    
    if [ \$MAF_LINES -le 1 ]; then
        # Empty MAF file
        touch ${prefix}.final.maf
    else
        Rscript /opt/annotate-maf-wrapper.R \\
            --hisens-rdata ${hisens_rdata} \\
            --maf-file ${maf} \\
            --output-dir annotation_output \\
            ${args}

        Rscript /opt/annotate-with-zygosity-germline.R \\
            --maf-file annotation_output/annotated.maf \\
            --output-maf ${prefix}.final.maf \\
            ${args}
    fi
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p annotation_output
    touch ${prefix}.final.maf
    """
}
