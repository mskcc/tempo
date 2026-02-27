process AGGREGATE_SOMATIC_HRDETECT {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(hrdetect_files)

    output:
    path("hrdetect.tsv"), emit: aggregated_hrdetect

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' ${hrdetect_files} > hrdetect.tsv
    """

    stub:
    """
    touch hrdetect.tsv
    """
}
