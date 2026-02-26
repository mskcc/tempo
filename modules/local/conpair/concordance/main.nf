process CONPAIR_CONCORDANCE {
    tag "$meta.id"
    label 'process_single'

    container "cmopipeline/conpair:v0.3.3"

    input:
    tuple val(meta), path(tumor_pileup), path(normal_pileup)

    output:
    tuple val(meta), path("*.concordance.tsv"), emit: concordance
    tuple val(meta), path("*.contamination.tsv"), emit: contamination
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    python /usr/local/bin/verify_concordance.py \\
        -T ${tumor_pileup} \\
        -N ${normal_pileup} \\
        -O ${prefix}.concordance.tsv

    python /usr/local/bin/estimate_tumor_normal_contamination.py \\
        -T ${tumor_pileup} \\
        -N ${normal_pileup} \\
        -O ${prefix}.contamination.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "0.3.3"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.concordance.tsv
    touch ${prefix}.contamination.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "stub"
    END_VERSIONS
    """
}
