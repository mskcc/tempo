process ASCAT_ALLELECOUNT {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/wtsicgp/ascatNgs:4.4.0'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path snp_gc_corrections

    output:
    tuple val(meta), path("ascat_alleleCount_*.tar.gz"), emit: alleles

    when:
    params.assay_type == 'genome'

    stub:
    """
    touch ascat_alleleCount_0.tar.gz
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p ascatResults
    export TMPDIR=\$(pwd)/tmp
    mkdir -p \$TMPDIR

    ascat.pl \\
        -o ./ascatResults \\
        -t ${tumor_bam} -n ${normal_bam} \\
        -sg ${snp_gc_corrections} \\
        -r ${fasta} \\
        -q 20 -g L \\
        -rs "homo_sapiens" -ra "GRCh37" -pr "WGS" \\
        -c ${task.cpus} \\
        -force \\
        -p allele_count

    tar -czf ascat_alleleCount_0.tar.gz ascatResults/
    """
}
