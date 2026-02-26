process SVABA_GERMLINE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/svaba:0.0.1' :
        'cmopipeline/svaba:0.0.1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("${prefix}.reheader.svaba.germline.sv.vcf.gz"), path("${prefix}.reheader.svaba.germline.sv.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def target_param = meta.target ? "-k ${meta.target}" : ""
    """
    svaba run \\
        -t ${bam} \\
        -G ${fasta} \\
        -p ${task.cpus * 2} \\
        -I \\
        -L 6 \\
        --id-string ${prefix} \\
        ${target_param} \\
        -z

    echo -e "${bam} ${meta.sample ?: meta.id}" > svaba.samplenames.tsv
    bcftools reheader \\
        --samples svaba.samplenames.tsv \\
        --output ${prefix}.reheader.svaba.germline.sv.vcf.gz \\
        ${prefix}.svaba.sv.vcf.gz
    bcftools index -f -t ${prefix}.reheader.svaba.germline.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: \$(svaba --version 2>&1 | head -1 | sed 's/.*version //')
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.reheader.svaba.germline.sv.vcf.gz
    touch ${prefix}.reheader.svaba.germline.sv.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svaba: 1.1.0
        bcftools: 1.20
    END_VERSIONS
    """
}
