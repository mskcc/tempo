process AGGREGATE_SOMATIC_SVSIGNATURES {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(svsignatures_files)

    output:
    path("sv_catalogues.pdf"), emit: catalogues
    path("sv_exposures.tsv"), emit: exposures

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Merge SV signature files (only sv_signatures sent by workflow)
    awk 'FNR==1 && NR!=1 {next} {print}' ${svsignatures_files} > sv_exposures.tsv

    # Create stub file for catalogues
    touch sv_catalogues.pdf
    """

    stub:
    """
    touch sv_catalogues.pdf sv_exposures.tsv
    """
}
