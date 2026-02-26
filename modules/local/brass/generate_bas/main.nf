process BRASS_GENERATE_BAS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/wtsicgp/pcap-core:5.5.0'

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.bas"), emit: bas

    when:
    params.assay_type == 'genome'

    stub:
    """
    touch ${meta.id}.bas
    """

    script:
    """
    bam_stats \\
        -i ${bam} \\
        -o ${meta.id}.bas \\
        -r ${fasta_fai} \\
        -@ ${task.cpus}
    """
}
