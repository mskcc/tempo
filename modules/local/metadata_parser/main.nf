process METADATA_PARSER {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/metadataparser:0.5.9' :
        'docker.io/cmopipeline/metadataparser:0.5.9' }"

    input:
    tuple val(meta), path(purity_out), path(maf_file), path(qc_output), path(msi_file), path(mutsig), path(polysolver_file)
    path(coding_bed)

    output:
    tuple val(meta), path("*.sample_data.txt"), emit: metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_metadata_file.py \\
        --sampleID ${prefix} \\
        --tumorID ${meta.tumor_id} \\
        --normalID ${meta.normal_id} \\
        --facetsPurity_out ${purity_out} \\
        --facetsQC ${qc_output} \\
        --MSIsensor_output ${msi_file} \\
        --mutational_signatures_output ${mutsig} \\
        --polysolver_output ${polysolver_file} \\
        --MAF_input ${maf_file} \\
        --coding_baits_BED ${coding_bed}

    mv ${prefix}_metadata.txt ${prefix}.sample_data.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metadataparser: 0.5.9
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sample_data.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metadataparser: 0.5.9
    END_VERSIONS
    """
}
