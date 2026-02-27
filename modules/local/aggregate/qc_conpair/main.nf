// Aggregate QC Conpair — matches original Tempo QcConpairAggregate
// Proper header construction and both concordance + contamination parsing
process AGGREGATE_QC_CONPAIR {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(concordance_files)
    path(contamination_files)

    output:
    path("concordance_qc.txt"), emit: concordance_qc
    path("contamination_qc.txt"), emit: contamination_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if ls *.concordance.txt 1> /dev/null 2>&1; then
        echo -e "Pair\\tConcordance" > concordance_qc.txt
        grep -v "concordance" *.concordance.txt | sed 's/.concordance.txt:/\\t/' | cut -f1,3 | sort -k1,1 >> concordance_qc.txt
    fi
    if ls *.contamination.txt 1> /dev/null 2>&1; then
        echo -e "Pair\\tSample_Type\\tSample_ID\\tContamination" > contamination_qc.txt
        grep -v "Contamination" *.contamination.txt | sed 's/.contamination.txt:/\\t/' | sort -k1,1 >> contamination_qc.txt
    fi
    touch concordance_qc.txt contamination_qc.txt
    """

    stub:
    """
    touch concordance_qc.txt contamination_qc.txt
    """
}
