// ASCAT run — matches original Tempo runAscat
// Genome-aware: GRCh37→37/HUMAN, GRCh38→38/HUMAN
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

    script:
    def genome = params.genome ?: 'GRCh37'
    def species = (genome in ['GRCh37', 'smallGRCh37', 'GRCh38']) ? 'HUMAN' : genome
    def assembly = (genome in ['GRCh38']) ? '38' : '37'
    """
    export TMPDIR=\$(pwd)/tmp
    mkdir -p \$TMPDIR

    for i in ascat_alleleCount_*.tar.gz ; do
        tar -xzf \$i
    done

    ascat.pl \\
        -o ./ascatResults \\
        -t ${tumor_bam} -n ${normal_bam} \\
        -sg ${snp_gc_corrections} \\
        -r ${fasta} \\
        -q 20 -g L \\
        -rs "${species}" -ra "${assembly}" -pr "WGS" \\
        -c ${task.cpus} \\
        -force
    """

    stub:
    """
    mkdir -p ascatResults
    touch ascatResults/sample.copynumber.caveman.csv
    touch ascatResults/sample.samplestatistics.txt
    """
}
