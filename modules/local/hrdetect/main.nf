process HRDETECT {
    tag "${meta.id}"
    label 'process_high'
    container 'cmopipeline/signaturetoolslib:0.0.1'

    input:
    tuple val(meta), path(maf_file), path(cnv_file), path(sv_file)
    path(hrdetect_script)

    output:
    tuple val(meta), path("*.hrdetect.tsv"), emit: hrdetect_output

    when:
    params.assay_type == 'genome'

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.hrdetect.tsv
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    # Build input TSV
    cat > hrdetect_input.tsv << INPUT_EOF
    sample_id\tmutation_file\tcnv_file\tsv_file
    ${meta.id}\t${maf_file}\t${cnv_file}\t${sv_file}
    INPUT_EOF

    # Run HRDetect analysis
    Rscript ${hrdetect_script} \\
        -i hrdetect_input.tsv \\
        -o ${prefix}.hrdetect.tsv \\
        -n ${task.cpus}
    """
}
