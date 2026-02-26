process AGGREGATE_SOMATIC_NETMHC {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(netmhc_files)

    output:
    path("mut_somatic_neoantigens.txt"), emit: aggregated_netmhc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' ${netmhc_files} > mut_somatic_neoantigens.txt
    """

    stub:
    """
    touch mut_somatic_neoantigens.txt
    """
}
