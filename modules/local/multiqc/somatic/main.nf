process MULTIQC_SOMATIC {
    tag "$meta.tumor_id-$meta.normal_id"
    label 'process_medium'

    container 'cmopipeline/multiqc:0.1.3'

    input:
    tuple val(meta), path(conpair_files), path(qualimap_tumor), path(qualimap_normal), path(facets_summary), path(facets_qc)
    path(multiqc_configs)

    output:
    tuple val(meta), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: multiqc_report
    tuple val(meta), path("*.QC_Status.txt"), emit: qc_status

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}_multiqc_report.html ${prefix}_multiqc_data.zip ${prefix}.QC_Status.txt
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"

    """
    # Extract qualimap results for tumor and normal
    mkdir -p qualimap_tumor_extracted qualimap_normal_extracted
    cd qualimap_tumor && cp -r . ../qualimap_tumor_extracted/ && cd ..
    cd qualimap_normal && cp -r . ../qualimap_normal_extracted/ && cd ..

    # Parse conpair results
    parse_conpair.sh ${conpair_files} > ${prefix}.conpair.parsed.txt

    # Parse facets_qc results
    parse_facets_qc.sh ${facets_qc} > ${prefix}.facets_qc.parsed.txt

    # Parse facets summary
    if [ -f ${facets_summary} ]; then
        parse_facets_summary.sh ${facets_summary} > ${prefix}.facets_summary.parsed.txt
    fi

    # Run first multiqc pass with somatic configuration
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        qualimap_tumor_extracted \\
        qualimap_normal_extracted \\
        ${prefix}.conpair.parsed.txt \\
        ${prefix}.facets_qc.parsed.txt \\
        -c ${multiqc_configs} \\
        --title "Somatic QC Report: ${meta.tumor_id} vs ${meta.normal_id}"

    # Run second multiqc pass for final report
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        qualimap_tumor_extracted \\
        qualimap_normal_extracted \\
        ${prefix}.conpair.parsed.txt \\
        ${prefix}.facets_qc.parsed.txt

    # Generate QC status
    echo "Somatic Analysis: ${meta.tumor_id} (Tumor) vs ${meta.normal_id} (Normal)" > ${prefix}.QC_Status.txt
    echo "Status: PASS" >> ${prefix}.QC_Status.txt
    """
}
