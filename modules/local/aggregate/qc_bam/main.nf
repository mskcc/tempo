// Aggregate QC BAM — matches original Tempo QcBamAggregate
// Passes assay type and cpus to R script
process AGGREGATE_QC_BAM {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/cmopipeline/alfred:v0.1.17'

    input:
    val(cohort)
    path(alfred_files)
    path(hsmetrics_files)

    output:
    path("alignment_qc.txt"), emit: alignment_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def assayType = params.assay_type == 'exome' ? 'wes' : 'wgs'
    """
    Rscript --no-init-file /usr/bin/create-aggregate-qc-file.R -n ${task.cpus} -a ${assayType}
    """

    stub:
    """
    touch alignment_qc.txt
    """
}
