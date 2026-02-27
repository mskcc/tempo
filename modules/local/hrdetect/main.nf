process HRDETECT {
    tag "${meta.id}"
    label 'process_high'
    container 'docker.io/cmopipeline/signaturetoolslib:0.0.1'

    input:
    tuple val(meta), path(maf_file), path(cnv_file), path(sv_file)
    path(hrdetect_script)

    output:
    tuple val(meta), path("*.hrdetect.tsv"), emit: hrdetect_output

    when:
    params.assay_type == 'genome'

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.hrdetect.tsv
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def genome = params.genome ?: 'GRCh37'
    def genome_version = genome == 'GRCh38' ? 'hg38' : 'hg19'
    """
    echo -e "sample\\tsv\\tmutations\\tcnv" > ${prefix}.tsv
    echo -e "${prefix}\\t${sv_file}\\t${maf_file}\\t${cnv_file}" >> ${prefix}.tsv
    Rscript ${hrdetect_script} ${prefix}.tsv ${genome_version} ${task.cpus}
    """
}
