process AGGREGATE_SOMATIC_METADATA {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(metadata_files)

    output:
    path("sample_data.txt"), emit: aggregated_metadata

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' ${metadata_files} > sample_data.txt
    """

    stub:
    """
    touch sample_data.txt
    """
}
