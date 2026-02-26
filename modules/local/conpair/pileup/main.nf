process CONPAIR_PILEUP {
    tag "$meta.id"
    label 'process_low'

    container "cmopipeline/conpair:v0.3.3"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.pileup"), emit: pileup
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /usr/local/bin/run_gatk_pileup_for_sample.py \\
        -B ${bam} \\
        -O ${prefix}.pileup \\
        --reference ${fasta} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "0.3.3"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "stub"
    END_VERSIONS
    """
}
