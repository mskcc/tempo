process SVCIRCOS {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/biocircos:0.0.1' :
        'docker.io/cmopipeline/biocircos:0.0.1' }"

    input:
    tuple val(meta), path(bedpe), path(cnv)
    val(genome)

    output:
    tuple val(meta), path("${prefix}.circos.html"), emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def genome_ = genome == 'GRCh37' ? 'hg19' : (genome == 'GRCh38' ? 'hg38' : genome)
    """
    Rscript /opt/biocircos/run_biocircos.R \\
        -b ${bedpe} \\
        -c ${cnv} \\
        -s ${prefix} \\
        -g ${genome_}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biocircos: 0.0.1
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.circos.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biocircos: 0.0.1
    END_VERSIONS
    """
}
