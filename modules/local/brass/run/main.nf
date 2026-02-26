process BRASS_RUN {
    tag "${meta.id}"
    label 'process_high'
    container 'cmopipeline/brass:0.0.2'

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

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    mkdir -p brass
    touch brass/out.vcf.gz
    touch brass/out.vcf.gz.tbi
    touch ${prefix}.brass.annot.vcf.gz
    touch ${prefix}.brass.annot.vcf.gz.tbi
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    # Assemble input data from BRASS_INPUT and BRASS_COVER directory outputs
    mkdir -p brass/tmpBrass
    cp -r ${brass_input_data}/tmpBrass/* brass/tmpBrass/
    cp -r ${brass_cover_data}/tmpBrass/* brass/tmpBrass/

    # Run BRASS full process
    brass.pl \\
        -p full \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -d ${brass_ref_dir} \\
        -v ${vagrent_ref_dir} \\
        -g ${fasta} \\
        -a ${ascat_samplestatistics} \\
        -i brass/tmpBrass \\
        -o brass

    # Reheader and annotate output
    bcftools reheader \\
        -h <(echo -e "##tumor_sample=${meta.tumor_id}\\n##normal_sample=${meta.normal_id}") \\
        brass/*.vcf.gz | \\
        bgzip -c > ${prefix}.brass.annot.vcf.gz

    tabix -p vcf ${prefix}.brass.annot.vcf.gz
    """
}
