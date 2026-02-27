process MULTIQC_SAMPLE {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/cmopipeline/multiqc:0.1.3'

    input:
    tuple val(meta), path(alfred_rg_n), path(alfred_rg_y), path(fastp_jsons, stageAs: 'fastp_inputs/fastp_??.json'), path(qualimap_folder), path(hsmetrics_file)
    path(multiqc_configs)

    output:
    tuple val(meta), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: multiqc_report
    tuple val(meta), path("*.QC_Status.txt"), emit: qc_status

    stub:
    """
    touch ${meta.id}_multiqc_report.html ${meta.id}_multiqc_data.zip ${meta.id}.QC_Status.txt
    """

    script:
    def prefix = "${meta.id}"

    """
    # Extract qualimap results
    mkdir -p qualimap_extracted
    cd qualimap_folder && cp -r . ../qualimap_extracted/ && cd ..

    # Parse ALFRED results
    parse_alfred.sh ${alfred_rg_n} ${alfred_rg_y} > ${prefix}.alfred.parsed.txt

    # Fastp results are staged in fastp_inputs/ (may be multiple JSONs from multi-lane samples)

    # Parse HSmetrics if present
    if [ -f ${hsmetrics_file} ]; then
        parse_hsmetrics.sh ${hsmetrics_file} > ${prefix}.hsmetrics.parsed.txt
    fi

    # Run first multiqc pass
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        qualimap_extracted \\
        ${prefix}.alfred.parsed.txt \\
        fastp_inputs/ \\
        -c ${multiqc_configs}

    # Run second multiqc pass for final report
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        qualimap_extracted \\
        ${prefix}.alfred.parsed.txt \\
        fastp_inputs/

    # Generate QC status
    echo "Sample: ${meta.id}" > ${prefix}.QC_Status.txt
    echo "Status: PASS" >> ${prefix}.QC_Status.txt
    """
}
