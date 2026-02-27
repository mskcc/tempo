// Aggregate Somatic LOHHLA — matches original Tempo SomaticAggregateLOHHLA
// Merges both HLAlossPrediction and IntegerCPN files
process AGGREGATE_SOMATIC_LOHHLA {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(lohhla_prediction_files)
    path(lohhla_cpn_files)

    output:
    path("HLAlossPrediction_CI.txt"), emit: hla_loss
    path("DNA.IntegerCPN_CI.txt"), emit: integer_cpn

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tmp ; TMPDIR=./tmp
    mkdir lohhla ; mv *.txt lohhla/ 2>/dev/null || true
    awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*HLAlossPrediction_CI.txt > HLAlossPrediction_CI.txt 2>/dev/null || touch HLAlossPrediction_CI.txt
    awk 'FNR==1 && NR!=1{next;}{print}' lohhla/*DNA.IntegerCPN_CI.txt > DNA.IntegerCPN_CI.txt 2>/dev/null || touch DNA.IntegerCPN_CI.txt
    """

    stub:
    """
    touch HLAlossPrediction_CI.txt DNA.IntegerCPN_CI.txt
    """
}
