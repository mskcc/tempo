process AGGREGATE_SOMATIC_SV {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
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
