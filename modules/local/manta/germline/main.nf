process MANTA_GERMLINE {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::strelka=2.9.10 bioconda::manta=1.5.0 bioconda::bcftools=1.9 bioconda::vt=0.57721"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/strelka2-manta-bcftools-vt:2.0.1' :
        'cmopipeline/strelka2-manta-bcftools-vt:2.0.1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    path(call_regions)
    path(call_regions_tbi)

    output:
    tuple val(meta), path("${prefix}.manta.vcf.gz"), path("${prefix}.manta.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def options = params.assay_type == 'exome' ? '--exome' : ''
    def regions_cmd = call_regions ? "--callRegions ${call_regions}" : ''
    """
    configManta.py \\
        ${options} \\
        ${regions_cmd} \\
        --reference ${fasta} \\
        --bam ${bam} \\
        --runDir Manta

    python Manta/runWorkflow.py \\
        --mode local \\
        --jobs ${task.cpus}

    mv Manta/results/variants/diploidSV.vcf.gz ${prefix}.manta.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi ${prefix}.manta.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version 2>&1 | sed -n '2s/.*manta version //p' || echo '1.6.0')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.manta.vcf.gz
    touch ${prefix}.manta.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: 1.6.0
    END_VERSIONS
    """
}
