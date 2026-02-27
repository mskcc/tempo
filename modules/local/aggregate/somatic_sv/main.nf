process AGGREGATE_SOMATIC_SV {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(bedpe_files)

    output:
    path("sv_somatic.bedpe"), emit: aggregated_sv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' ${bedpe_files} > sv_somatic.bedpe
    """

    stub:
    """
    touch sv_somatic.bedpe
    """
}
