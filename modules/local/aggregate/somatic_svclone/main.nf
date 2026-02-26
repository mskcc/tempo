process AGGREGATE_SOMATIC_SVCLONE {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(svclone_files)

    output:
    path("svclone_sv_cluster_certainty.tsv"), emit: sv_clusters
    path("svclone_snv_cluster_certainty.tsv"), emit: snv_clusters

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Merge cluster certainty files (workflow collects both SV and SNV together)
    awk 'FNR==1 && NR!=1 {next} {print}' ${svclone_files} > svclone_sv_cluster_certainty.tsv

    # Create stub file for SNV clusters
    touch svclone_snv_cluster_certainty.tsv
    """

    stub:
    """
    touch svclone_sv_cluster_certainty.tsv svclone_snv_cluster_certainty.tsv
    """
}
