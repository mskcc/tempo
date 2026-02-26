process MUTSIG {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/temposig:0.2.3'

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.mutsig.txt"), emit: mutsig_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def cosmic_version = params.cosmic_version ?: '92'
    """
    mkdir -p mutsig_output

    # Convert MAF to catalog format using maf2cat2.R
    Rscript /opt/maf2cat2.R \\
        --maf-file ${maf} \\
        --output-dir mutsig_output \\
        ${args}

    # Run tempoSig analysis with cosmic version
    Rscript /opt/tempoSig.R \\
        --catalog-file mutsig_output/catalog.txt \\
        --cosmic-version ${cosmic_version} \\
        --output-file ${prefix}.mutsig.txt \\
        ${args}
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p mutsig_output
    touch ${prefix}.mutsig.txt
    """
}
