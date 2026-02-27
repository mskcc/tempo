// BRASS run — matches original Tempo SomaticRunBRASS
// Genome-aware: GRCh37→37/HUMAN, GRCh38→38/HUMAN
process BRASS_RUN {
    tag "${meta.id}"
    label 'process_high'
    container 'docker.io/cmopipeline/brass:0.0.2'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(tumor_bas), path(normal_bam), path(normal_bai), path(normal_bas), path(brass_input_data), path(brass_cover_data), path(ascat_samplestatistics)
    path(fasta)
    path(fasta_fai)
    path(brass_ref_dir)
    path(vagrent_ref_dir)

    output:
    tuple val(meta), path("brass/*.{vcf.gz,vcf.gz.tbi}"), emit: brass_output
    tuple val(meta), path("*.brass.annot.vcf.gz"), path("*.brass.annot.vcf.gz.tbi"), emit: brass_vcf

    when:
    params.assay_type == 'genome'

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def genome = params.genome ?: 'GRCh37'
    def species = (genome in ['GRCh37', 'smallGRCh37', 'GRCh38']) ? 'HUMAN' : genome
    def assembly = (genome in ['GRCh38']) ? '38' : '37'
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR

    # Assemble input data from BRASS_INPUT and BRASS_COVER
    mkdir -p brass/tmpBrass
    cp -r ${brass_input_data}/tmpBrass/* brass/tmpBrass/
    cp -r ${brass_cover_data}/tmpBrass/* brass/tmpBrass/

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
        -ss ${ascat_samplestatistics} \\
        -o brass \\
        -p full

    # Reheader with proper sample names
    echo -e "TUMOUR ${meta.tumor_id}\\nNORMAL ${meta.normal_id}" > samples.txt
    bcftools reheader \\
        --samples samples.txt \\
        brass/*.vcf.gz | \\
    bcftools sort | \\
    bgzip -c > ${prefix}.brass.annot.vcf.gz

    tabix -p vcf ${prefix}.brass.annot.vcf.gz
    """

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p brass
    touch brass/out.vcf.gz
    touch brass/out.vcf.gz.tbi
    touch ${prefix}.brass.annot.vcf.gz
    touch ${prefix}.brass.annot.vcf.gz.tbi
    """
}
