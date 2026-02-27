// Aggregate Somatic SVclone — matches original Tempo SomaticAggregateSVclone
// Merges both SV and SNV cluster certainty files
process AGGREGATE_SOMATIC_SVCLONE {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    val(cohort)
    path(sv_files, stageAs: 'sv/*')
    path(snv_files, stageAs: 'snv/*')

    output:
    path("svclone_sv_cluster_certainty.tsv"), emit: sv_clusters
    path("svclone_snv_cluster_certainty.tsv"), emit: snv_clusters

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' sv/* > svclone_sv_cluster_certainty.tsv
    awk 'FNR==1 && NR!=1{next;}{print}' snv/* > svclone_snv_cluster_certainty.tsv
    """

    stub:
    """
    touch svclone_sv_cluster_certainty.tsv svclone_snv_cluster_certainty.tsv
    """
}
