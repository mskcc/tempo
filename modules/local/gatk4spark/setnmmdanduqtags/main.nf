process GATK4SPARK_SETNMMDANDUQTAGS {
    tag "$meta.id"
    label 'process_low'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4-spark:4.6.1.0--hdfd78af_0':
        'biocontainers/gatk4-spark:4.6.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam, name:"input/*"), path(bai, name:"input/*"), path(intervals)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.bam")	 	 , emit: bam
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[gatk SetNmMdAndUqTags] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk \\
        --java-options "-Dsamjdk.compression_level=2 -Xmx${avail_mem}M -XX:-UsePerfData" \\
        PrintReadsSpark \\
        --input ${bam} \\
        --intervals ${intervals} \\
        --reference ${fasta} \\
        --create-output-bam-index true \\
        --spark-master local[${task.cpus}] \\
        --output tmp.bam
    gatk \\
        --java-options "-Dsamjdk.compression_level=2 -Dsamjdk.compression_level=1 -Xmx${avail_mem}M" \\
        SetNmMdAndUqTags \\
        $args \\
        --CREATE_INDEX false \\
        --INPUT tmp.bam \\
        --OUTPUT ${prefix}.bam \\
        --REFERENCE_SEQUENCE ${fasta}
    rm -rf tmp.bam tmp.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$( echo \$(gatk SetNmMdAndUqTags --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$( echo \$(gatk SetNmMdAndUqTags --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
