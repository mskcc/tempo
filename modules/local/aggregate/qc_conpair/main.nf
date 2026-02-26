process AGGREGATE_QC_CONPAIR {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(concordance_files)

    output:
    path("concordance_qc.txt"), emit: concordance_qc
    path("contamination_qc.txt"), emit: contamination_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Parse concordance files into aggregated TSV (only concordance sent by workflow)
    awk 'FNR==1 && NR!=1 {next} {print}' ${concordance_files} > concordance_qc.txt

    # Create stub file for contamination
    touch contamination_qc.txt
    """

    stub:
    """
    touch concordance_qc.txt contamination_qc.txt
    """
}
