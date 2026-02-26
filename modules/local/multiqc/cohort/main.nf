process MULTIQC_COHORT {
    tag "cohort"
    label 'process_low'

    container 'cmopipeline/multiqc:0.1.3'

    input:
    path(qc_files)

    output:
    path("*multiqc_report*.html"), emit: report
    path("*multiqc_data*"),        emit: data
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "cohort"
    """
    multiqc \\
        -n ${prefix}_multiqc_report \\
        -o . \\
        --title "Cohort QC Report" \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/.*version //')
    END_VERSIONS
    """

    stub:
    """
    touch cohort_multiqc_report.html cohort_multiqc_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.14
    END_VERSIONS
    """
}
