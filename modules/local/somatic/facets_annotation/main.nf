process SOMATIC_FACETS_ANNOTATION {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/facets-suite-preview-htstools:0.0.1'

    input:
    tuple val(meta), path(hisens_rdata), path(maf)

    output:
    tuple val(meta), path("*.somatic.final.maf"), emit: final_maf
    path("file-size.txt"), emit: maf_size

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
        touch ${prefix}.somatic.final.maf
        echo "0" > file-size.txt
    else
        Rscript /opt/annotate-maf-wrapper.R \\
            --hisens-rdata ${hisens_rdata} \\
            --maf-file ${maf} \\
            --output-dir annotation_output \\
            ${args}

        Rscript /opt/annotate-with-zygosity-somatic.R \\
            --maf-file annotation_output/annotated.maf \\
            --output-maf ${prefix}.somatic.final.maf \\
            ${args}

        wc -l < ${prefix}.somatic.final.maf > file-size.txt
    fi
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p annotation_output
    touch ${prefix}.somatic.final.maf file-size.txt
    """
}
