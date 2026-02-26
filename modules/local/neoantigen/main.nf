process NEOANTIGEN {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/neoantigen:0.3.3'

    input:
    tuple val(meta), path(polysolver_file), path(maf_file)
    path(neoantigen_cdna)
    path(neoantigen_cds)

    output:
    tuple val(meta), path("*.all_neoantigen_predictions.txt"), emit: predictions
    tuple val(meta), path("*neoantigens.maf"), emit: neoantigen_maf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p neoantigen_output

    python /opt/neoantigen.py \\
        --polysolver-file ${polysolver_file} \\
        --maf-file ${maf_file} \\
        --cdna-fasta ${neoantigen_cdna} \\
        --cds-fasta ${neoantigen_cds} \\
        --output-dir neoantigen_output \\
        ${args}

    # Prepend sample ID to predictions
    awk -v sample="${meta.tumor_id}" 'NR==1 {print "sample_id," \$0} NR>1 {print sample "," \$0}' \\
        neoantigen_output/predictions.txt > ${prefix}.all_neoantigen_predictions.txt

    # Create neoantigen MAF file with sample ID
    awk -v sample="${meta.tumor_id}" 'NR==1 {print} NR>1 {print}' \\
        neoantigen_output/neoantigens.maf > ${prefix}.neoantigens.maf
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p neoantigen_output
    touch ${prefix}.all_neoantigen_predictions.txt ${prefix}.neoantigens.maf
    """
}
