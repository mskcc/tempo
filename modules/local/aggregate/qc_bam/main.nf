process AGGREGATE_QC_BAM {
    label 'process_single'
    container 'cmopipeline/alfred:v0.1.17'

    input:
    path(alfred_files)
    path(hsmetrics_files)

    output:
    path("alignment_qc.txt"), emit: alignment_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript create-aggregate-qc-file.R
    """

    stub:
    """
    touch alignment_qc.txt
    """
}
