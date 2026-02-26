process AGGREGATE_SOMATIC_LOHHLA {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(lohhla_files)

    output:
    path("HLAlossPrediction_CI.txt"), emit: hla_loss
    path("DNA.IntegerCPN_CI.txt"), emit: integer_cpn

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Merge HLA loss prediction files (only predictions sent by workflow)
    awk 'FNR==1 && NR!=1 {next} {print}' ${lohhla_files} > HLAlossPrediction_CI.txt

    # Create stub file for integer CPN
    touch DNA.IntegerCPN_CI.txt
    """

    stub:
    """
    touch HLAlossPrediction_CI.txt DNA.IntegerCPN_CI.txt
    """
}
