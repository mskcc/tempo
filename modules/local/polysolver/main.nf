process POLYSOLVER {
    tag "$meta.id"
    label 'process_high'

    container "sachet/polysolver:v4"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.hla.txt"),    emit: hla_types
    tuple val(meta), path("*"),            emit: results
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /home/polysolver/scripts/shell_call_hla_type \\
        ${bam} \\
        Unknown \\
        1 \\
        hg19 \\
        STRELKA \\
        0 \\
        ${prefix}

    mv ${prefix}/winners.hla.txt ${prefix}.hla.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: "v4"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.hla.txt
    touch ${prefix}_winners.hla.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polysolver: "stub"
    END_VERSIONS
    """
}
