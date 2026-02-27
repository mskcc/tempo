// Aggregate Somatic MAF — matches original Tempo SomaticAggregateMaf
process AGGREGATE_SOMATIC_MAF {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(maf_files)

    output:
    path("mut_somatic.maf"), emit: aggregated_maf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tmp ; TMPDIR=./tmp
    mkdir mut ; mv *.maf mut/
    for i in mut/*.maf ; do
        if [ \$( cat \$i | wc -l ) -gt 1 ] ; then cat \$i ; fi
    done | grep ^Hugo_Symbol | head -n 1 > mut_somatic.maf
    cat mut/*.maf | grep -Ev "^#|^Hugo_Symbol" | sort -k5,5V -k6,6n >> mut_somatic.maf
    """

    stub:
    """
    touch mut_somatic.maf
    """
}
