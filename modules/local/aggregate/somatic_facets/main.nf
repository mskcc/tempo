process AGGREGATE_SOMATIC_FACETS {
    tag "${cohort}"
    label 'process_single'
    container 'docker.io/library/ubuntu:22.04'

    input:
    tuple val(cohort), path(purity), path(Hisens), path(outLog), path(armLev), path(geneLev)

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
    mkdir facets_tmp
    mv *_OUT.txt facets_tmp/
    mv *{purity,hisens}.seg facets_tmp/
    awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_hisens.seg > cna_hisens_run_segmentation.seg
    awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_purity.seg > cna_purity_run_segmentation.seg
    awk 'FNR==1 && NR!=1{next;}{print}' facets_tmp/*_OUT.txt > cna_facets_run_info.txt
    mv *{gene_level,arm_level}.txt facets_tmp/
    cat facets_tmp/*gene_level.txt | head -n 1 > cna_genelevel.txt
    awk -v FS='\t' '{ if (\$24 != "DIPLOID" && (\$25 == "PASS" || \$25 == "RESCUE" ))  print \$0 }' facets_tmp/*gene_level.txt >> cna_genelevel.txt
    cat facets_tmp/*arm_level.txt | head -n 1 > cna_armlevel.txt
    cat facets_tmp/*arm_level.txt | grep -v "DIPLOID" | grep -v "Tumor_Sample_Barcode" >> cna_armlevel.txt || [[ \$? == 1 ]]
    """

    stub:
    """
    touch cna_hisens_run_segmentation.seg cna_purity_run_segmentation.seg cna_armlevel.txt cna_genelevel.txt cna_facets_run_info.txt
    """
}
