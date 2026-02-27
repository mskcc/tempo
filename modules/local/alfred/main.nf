process ALFRED {
    tag "${meta.id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/alfred:v0.1.17'

    input:
    tuple val(meta), path(bam), path(bai), path(targets), path(targets_index)
    path fasta

    output:
    tuple val(meta), path("*.alfred.per_readgroup.tsv.gz"), path("*.alfred.tsv.gz"), emit: alfred_qc
    path "versions.yml", emit: versions

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}.alfred.tsv.gz
    touch ${prefix}.alfred.per_readgroup.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alfred: 0.1.17
    END_VERSIONS
    """

    script:
    def prefix = "${meta.id}"
    def targets_arg = targets && targets.name != 'NO_FILE' ? "--bed ${targets}" : ""
    """
    alfred qc \\
        --reference ${fasta} \\
        ${targets_arg} \\
        --outfile ${prefix}.alfred.tsv.gz \\
        ${bam} && \\
    Rscript --no-init-file /opt/alfred/scripts/stats.R ${prefix}.alfred.tsv.gz

    alfred qc \\
        --reference ${fasta} \\
        ${targets_arg} \\
        --ignore \\
        --outfile ${prefix}.alfred.per_readgroup.tsv.gz \\
        ${bam} && \\
    Rscript --no-init-file /opt/alfred/scripts/stats.R ${prefix}.alfred.per_readgroup.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alfred: \$(alfred version 2>&1 | head -1 | sed 's/.*v//')
    END_VERSIONS
    """
}
