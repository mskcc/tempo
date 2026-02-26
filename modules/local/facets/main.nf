process FACETS {
    tag "$meta.id"
    label 'process_medium'

    container "cmopipeline/facets-suite-preview-htstools:0.0.1"

    input:
    tuple val(meta), path(snp_pileup)

    output:
    tuple val(meta), path("*_hisens.Rdata"),                   emit: hisens_rdata
    tuple val(meta), path("*_hisens.seg"),                     emit: hisens_seg
    tuple val(meta), path("*_purity.Rdata"),                   emit: purity_rdata, optional: true
    tuple val(meta), path("*_purity.seg"),                     emit: purity_seg, optional: true
    tuple val(meta), path("*.out"),                            emit: purity
    tuple val(meta), path("*{.png,.pdf}"),                     emit: plots,     optional: true
    tuple val(meta), path("*.facets_qc.txt"),                  emit: facets_qc, optional: true
    tuple val(meta), path("*.arm_level.txt"),                  emit: arm_level, optional: true
    tuple val(meta), path("*.gene_level.txt"),                 emit: gene_level, optional: true
    tuple val(meta), path("*.txt"),                            emit: summary,   optional: true
    tuple val(meta), path("*"),                                emit: facets_output
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    touch .Rprofile

    Rscript /usr/bin/facets-suite/run-facets-wrapper.R \\
        --counts-file ${snp_pileup} \\
        --sample-id ${prefix} \\
        --directory . \\
        --everything \\
        --legacy-output T \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: \$(Rscript -e 'cat(as.character(packageVersion("facets")))' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}_hisens.Rdata
    touch ${prefix}_hisens.seg
    touch ${prefix}_purity.Rdata
    touch ${prefix}_purity.seg
    touch ${prefix}.out
    touch ${prefix}.png
    touch ${prefix}.facets_qc.txt
    touch ${prefix}.arm_level.txt
    touch ${prefix}.gene_level.txt
    touch ${prefix}.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: "stub"
    END_VERSIONS
    """
}
