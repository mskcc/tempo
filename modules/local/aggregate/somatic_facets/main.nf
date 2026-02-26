process AGGREGATE_SOMATIC_FACETS {
    label 'process_single'
    container 'ubuntu:22.04'

    input:
    path(facets_files)

    output:
    path("cna_hisens_run_segmentation.seg"), emit: hisens_seg
    path("cna_purity_run_segmentation.seg"), emit: purity_seg
    path("cna_armlevel.txt"), emit: arm_level
    path("cna_genelevel.txt"), emit: gene_level
    path("cna_facets_run_info.txt"), emit: facets_info

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Merge hisens segmentation files (only hisens_seg sent by workflow)
    awk 'FNR==1 && NR!=1 {next} {print}' ${facets_files} > cna_hisens_run_segmentation.seg

    # Create stub files for other outputs
    touch cna_purity_run_segmentation.seg cna_armlevel.txt cna_genelevel.txt cna_facets_run_info.txt
    """

    stub:
    """
    touch cna_hisens_run_segmentation.seg cna_purity_run_segmentation.seg cna_armlevel.txt cna_genelevel.txt cna_facets_run_info.txt
    """
}
