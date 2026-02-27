// BRASS input — matches original Tempo SomaticRunBRASSInput
// Genome-aware: GRCh37→37/HUMAN, GRCh38→38/HUMAN
process BRASS_INPUT {
    tag "${meta.id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/brass:0.0.2'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(tumor_bas), path(normal_bam), path(normal_bai), path(normal_bas)
    path fasta
    path fasta_fai
    path brass_ref_dir
    path vagrent_ref_dir

    output:
    tuple val(meta), path("brass_input_data"), emit: brass_input

    when:
    params.assay_type == 'genome'

    script:
    def genome = params.genome ?: 'GRCh37'
    def species = (genome in ['GRCh37', 'smallGRCh37', 'GRCh38']) ? 'HUMAN' : genome
    def assembly = (genome in ['GRCh38']) ? '38' : '37'
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass
    for i in rho Ploidy GenderChr GenderChrFound ; do echo \$i ; done > samplestatistics.txt

    brass.pl -j 4 -k 4 -c ${task.cpus} \\
        -d ${brass_ref_dir}/HiDepth.bed.gz \\
        -f ${brass_ref_dir}/brass_np.groups.gz \\
        -g ${fasta} \\
        -s "${species}" -as "${assembly}" -pr "WGS" \\
        -g_cache ${vagrent_ref_dir}/vagrent.cache.gz \\
        -vi ${brass_ref_dir}/viral.genomic.fa.2bit \\
        -mi ${brass_ref_dir}/all_ncbi_bacteria \\
        -b ${brass_ref_dir}/500bp_windows.gc.bed.gz \\
        -ct ${brass_ref_dir}/CentTelo.tsv \\
        -cb ${brass_ref_dir}/cytoband.txt \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -ss samplestatistics.txt \\
        -o brass \\
        -p input

    mv brass brass_input_data
    """

    stub:
    """
    mkdir -p brass_input_data/tmpBrass/progress
    touch brass_input_data/tmpBrass/input.0
    touch brass_input_data/tmpBrass/progress/input.0.done
    """
}
