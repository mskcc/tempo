process GATK4_MERGEMUTECTSTATS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data' :
        'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3' }"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${prefix}.merged.stats"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = stats.collect { "--stats ${it}" }.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MergeMutectStats] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MergeMutectStats \\
        ${input_list} \\
        --output ${prefix}.merged.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n '/GATK.*v/s/.*v//p')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.1.0.0
    END_VERSIONS
    """
}
