process MULTIQC_COHORT {
    tag "$cohort"
    label 'process_low'

    container 'cmopipeline/multiqc:0.1.3'

    input:
    tuple val(cohort), path(fastp_tumor), path(fastp_normal), path(alfred_ign_y_tumor), path(alfred_ign_y_normal), path(alfred_ign_n_tumor), path(alfred_ign_n_normal), path(concord_file), path(contami_file), path(facets_summary), path(facets_qc), path(qualimap_tumor), path(qualimap_normal), path(hsmetrics_tumor), path(hsmetrics_normal), path(multiqc_configs)

    output:
    tuple val(cohort), path("*multiqc_report*.html"), path("*multiqc_data*.zip"), emit: cohort_multiqc_report

    stub:
    """
    touch ${cohort}_multiqc_report.html ${cohort}_multiqc_data.zip
    """

    script:
    def prefix = "${cohort}"

    """
    # Create directory structure for all inputs
    mkdir -p cohort_data/fastp cohort_data/alfred cohort_data/conpair cohort_data/facets cohort_data/qualimap cohort_data/hsmetrics

    # Organize fastp results
    cp ${fastp_tumor} cohort_data/fastp/tumor_fastp.json
    cp ${fastp_normal} cohort_data/fastp/normal_fastp.json

    # Organize ALFRED results
    cp ${alfred_ign_y_tumor} cohort_data/alfred/tumor_ign_y.tsv.gz
    cp ${alfred_ign_y_normal} cohort_data/alfred/normal_ign_y.tsv.gz
    cp ${alfred_ign_n_tumor} cohort_data/alfred/tumor_ign_n.tsv.gz
    cp ${alfred_ign_n_normal} cohort_data/alfred/normal_ign_n.tsv.gz

    # Organize conpair results
    cp ${concord_file} cohort_data/conpair/concordance.txt
    cp ${contami_file} cohort_data/conpair/contamination.txt

    # Organize facets results
    cp ${facets_summary} cohort_data/facets/summary.txt
    cp ${facets_qc} cohort_data/facets/qc.txt

    # Organize qualimap results
    mkdir -p cohort_data/qualimap/tumor cohort_data/qualimap/normal
    cd qualimap_tumor && cp -r . ../cohort_data/qualimap/tumor/ && cd ..
    cd qualimap_normal && cp -r . ../cohort_data/qualimap/normal/ && cd ..

    # Organize HSmetrics results
    cp ${hsmetrics_tumor} cohort_data/hsmetrics/tumor_metrics.txt
    cp ${hsmetrics_normal} cohort_data/hsmetrics/normal_metrics.txt

    # Parse all data for multiqc
    parse_alfred_cohort.sh cohort_data/alfred > ${prefix}.alfred.parsed.txt
    parse_conpair_cohort.sh cohort_data/conpair > ${prefix}.conpair.parsed.txt
    parse_facets_cohort.sh cohort_data/facets > ${prefix}.facets.parsed.txt
    parse_hsmetrics_cohort.sh cohort_data/hsmetrics > ${prefix}.hsmetrics.parsed.txt

    # Run multiqc for complete cohort report
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        cohort_data \\
        ${prefix}.alfred.parsed.txt \\
        ${prefix}.conpair.parsed.txt \\
        ${prefix}.facets.parsed.txt \\
        ${prefix}.hsmetrics.parsed.txt \\
        -c ${multiqc_configs} \\
        --title "Cohort QC Report: ${cohort}"

    # Run final multiqc pass
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        cohort_data \\
        ${prefix}.alfred.parsed.txt \\
        ${prefix}.conpair.parsed.txt \\
        ${prefix}.facets.parsed.txt \\
        ${prefix}.hsmetrics.parsed.txt
    """
}
