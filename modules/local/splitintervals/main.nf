process SPLIT_INTERVALS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://broadinstitute/gatk:4.1.0.0' :
        'broadinstitute/gatk:4.1.0.0' }"

    input:
    tuple val(meta), path(fasta), path(fai), path(dict)
    path(intervals)
    val(scatter_count)

    output:
    path("*.interval_list"), emit: interval_lists
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.assay_type == 'genome' ? 'INTERVAL_SUBDIVISION' : 'BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW'
    """
    gatk SplitIntervals \\
        --reference ${fasta} \\
        --intervals ${intervals} \\
        --scatter-count ${scatter_count} \\
        --subdivision-mode ${mode} \\
        --output scattered

    for i in scattered/*.interval_list; do
        BASENAME=\$(basename \$i)
        mv \$i scattered-\$BASENAME
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '2s/.*v//p')
    END_VERSIONS
    """

    stub:
    """
    touch scattered-0001.interval_list scattered-0002.interval_list
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.1.0.0
    END_VERSIONS
    """
}
