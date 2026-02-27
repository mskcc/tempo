// FACETS — matches original Tempo DoFacets
// Includes retry logic (4 attempts with seed variation), post-processing
process FACETS {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/cmopipeline/facets-suite-preview-htstools:0.0.1"

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
    tuple val(meta), path("*_OUT.txt"),                        emit: summary_out, optional: true
    tuple val(meta), path("*_hisens_samplestatistics.txt"),    emit: hisens_stats, optional: true
    tuple val(meta), path("*_purity_samplestatistics.txt"),    emit: purity_stats, optional: true
    tuple val(meta), path("*"),                                emit: facets_output
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}__${meta.normal_id}"
    def facets_cval = params.facets_cval ?: '100'
    def facets_snp_nbhd = params.facets_snp_nbhd ?: '250'
    def facets_ndepth = params.facets_ndepth ?: '35'
    def facets_min_nhet = params.facets_min_nhet ?: '25'
    def facets_purity_cval = params.facets_purity_cval ?: '500'
    def facets_purity_min_nhet = params.facets_purity_min_nhet ?: '25'
    def facets_genome = params.facets_genome ?: 'hg19'
    def facets_seed = params.facets_seed ?: '100'
    """
    touch .Rprofile

    mkdir -p facets_output

    # Retry logic: up to 4 attempts with seed variation (matches original Tempo)
    set +e
    i=1
    seed=\$((${facets_seed}-1))
    attemptNumber=0

    while [ \$i -eq 1 ]
    do
        attemptNumber=\$(( attemptNumber + 1 ))
        if [ \$attemptNumber -gt 4 ]; then
            break
        fi
        seed=\$((seed+i))

        Rscript /usr/bin/facets-suite/run-facets-wrapper.R \\
            --cval ${facets_cval} \\
            --snp-window-size ${facets_snp_nbhd} \\
            --normal-depth ${facets_ndepth} \\
            --min-nhet ${facets_min_nhet} \\
            --purity-cval ${facets_purity_cval} \\
            --purity-min-nhet ${facets_purity_min_nhet} \\
            --genome ${facets_genome} \\
            --counts-file ${snp_pileup} \\
            --sample-id ${prefix} \\
            --directory facets_output \\
            --facets-lib-path /usr/local/lib/R/site-library \\
            --seed \$seed \\
            --everything \\
            --legacy-output T \\
            ${args}

        i=\$?
    done
    set -e

    # Post-processing: summarize project (matches original Tempo)
    python3 /usr/bin/summarize_project.py \\
        -p ${prefix} \\
        -c facets_output/*cncf.txt \\
        -o facets_output/*out \\
        -s facets_output/*seg || true

    # Generate sample statistics for hisens and purity (matches original Tempo)
    if [ -f facets_output/${prefix}_hisens.Rdata ] ; then
        Rscript /usr/bin/facets-suite/generate_samplestatistics.R \\
            facets_output/${prefix}_hisens.Rdata \\
            ${prefix}_hisens || true
    fi

    if [ -f facets_output/${prefix}_purity.Rdata ] ; then
        Rscript /usr/bin/facets-suite/generate_samplestatistics.R \\
            facets_output/${prefix}_purity.Rdata \\
            ${prefix}_purity || true
    fi

    # Move all outputs to working directory
    mv facets_output/* . 2>/dev/null || true

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
    touch ${prefix}_OUT.txt
    touch ${prefix}_hisens_samplestatistics.txt
    touch ${prefix}_purity_samplestatistics.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        facets: "stub"
    END_VERSIONS
    """
}
