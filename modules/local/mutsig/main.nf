// MUTSIG / TempoSig — matches original Tempo RunMutationSignatures
// Two-stage: maf2cat2.R → tempoSig.R with fixed statistical params
process MUTSIG {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/temposig:0.2.3'

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.mutsig.txt"), emit: mutsig_results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def cosmic = params.cosmic ?: 'v3'
    """
    maf2cat2.R ${maf} \\
        ${prefix}.trinucmat.txt

    tempoSig.R --cosmic_${cosmic} --pvalue --nperm 10000 --seed 132 ${prefix}.trinucmat.txt \\
        ${prefix}.mutsig.txt
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.mutsig.txt
    """
}
