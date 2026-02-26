process BRASS_COVER {
    tag "${meta.id}"
    label 'process_medium'
    container 'cmopipeline/brass:0.0.2'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(tumor_bas), path(normal_bam), path(normal_bai), path(normal_bas)
    path fasta
    path fasta_fai
    path brass_ref_dir
    path vagrent_ref_dir

    output:
    tuple val(meta), path("brass_cover_data"), emit: brass_cover

    when:
    params.assay_type == 'genome'

    stub:
    """
    mkdir -p brass_cover_data/tmpBrass/cover brass_cover_data/tmpBrass/progress
    touch brass_cover_data/tmpBrass/cover/cover.0
    touch brass_cover_data/tmpBrass/progress/cover.0.done
    """

    script:
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass
    brass.pl \\
        -p cover \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -d ${brass_ref_dir} \\
        -v ${vagrent_ref_dir} \\
        -g ${fasta} \\
        -o brass

    # Move output to uniquely-named directory to avoid staging collisions
    mv brass brass_cover_data
    """
}
