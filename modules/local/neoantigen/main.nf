// Neoantigen prediction — matches original Tempo RunNeoantigen
// Includes TMPDIR setup, config file, tab-separated output
process NEOANTIGEN {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/neoantigen:0.3.3'

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
    def output_dir = "neoantigen"
    """
    export TMPDIR=\$PWD/${output_dir}-tmp/
    mkdir -p ${output_dir}-tmp
    chmod 777 ${output_dir}-tmp

    python /usr/local/bin/neoantigen/neoantigen.py \\
        --config_file /usr/local/bin/neoantigen/neoantigen-docker.config \\
        --sample_id ${prefix} \\
        --hla_file ${polysolver_file} \\
        --maf_file ${maf_file} \\
        --threads ${task.cpus} \\
        --output_dir ${output_dir} \\
        ${args}

    awk 'NR==1 {printf("%s\\t%s\\n", "sample", \$0)} NR>1 {printf("%s\\t%s\\n", "${prefix}", \$0) }' ${output_dir}/*.all_neoantigen_predictions.txt > ${prefix}.all_neoantigen_predictions.txt

    cp ${output_dir}/${prefix}.neoantigens.maf ${prefix}.neoantigens.maf || true
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p neoantigen
    touch ${prefix}.all_neoantigen_predictions.txt ${prefix}.neoantigens.maf
    """
}
