process ASCAT_RUN {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/wtsicgp/ascatNgs:4.4.0'

    input:
    tuple val(meta), path(ascat_tars), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(fasta)
    path(fasta_fai)
    path(snp_gc_corrections)

    output:
    tuple val(meta), path("ascatResults/*.copynumber.caveman.csv"), emit: caveman
    tuple val(meta), path("ascatResults/*.samplestatistics.txt"), emit: ascat_samplestatistics

    when:
    params.assay_type == 'genome'

    stub:
    """
    mkdir -p ascatResults
    touch ascatResults/sample.copynumber.caveman.csv
    touch ascatResults/sample.samplestatistics.txt
    """

    script:
    """
    # Extract all tar files
    mkdir -p ascat_work
    for tar in ${ascat_tars}; do
        tar -xzf \$tar -C ascat_work
    done

    # Run ASCAT full process
    ascat.pl \\
        full \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -r ${fasta} \\
        -g ${snp_gc_corrections} \\
        -i ascat_work \\
        -o ascatResults
    """
}
