process SV_SIGNATURES {
    tag "${meta.id}"
    label 'process_high'
    container 'cmopipeline/signaturetoolslib:0.0.1'

    input:
    tuple val(meta), path(bedpe)
    path(sv_signature_script)

    output:
    tuple val(meta), path("*_catalogues.pdf"), emit: catalogues
    tuple val(meta), path("*_exposures.tsv"), emit: sv_signatures

    when:
    params.assay_type == 'genome'

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}_catalogues.pdf
    touch ${prefix}_exposures.tsv
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def genome_version = params.genome == 'GRCh38' ? 'hg38' : 'hg19'
    """
    Rscript ${sv_signature_script} \\
        -i ${bedpe} \\
        -g ${genome_version} \\
        -n ${task.cpus} \\
        -s ${prefix} \\
        -o .
    """
}
