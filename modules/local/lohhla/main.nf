process LOHHLA {
    tag "$meta.id"
    label 'process_high'

    container "cmopipeline/lohhla:1.1.7"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(hla_types)
    path hla_fasta
    path hla_dat

    output:
    tuple val(meta), path("*.DNA.HLAlossPrediction_CI.xls"), emit: predictions, optional: true
    tuple val(meta), path("*"),                               emit: results
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript /opt/lohhla/LOHHLAscript.R \\
        --patientId ${prefix} \\
        --normalBAMfile ${normal_bam} \\
        --tumorBAMfile ${tumor_bam} \\
        --HLAfastaLoc ${hla_fasta} \\
        --HLAexonLoc ${hla_dat} \\
        --hlaPath . \\
        --CopyNumLoc . \\
        --mappingStep TRUE \\
        --minCoverageFilter 10 \\
        --fishingStep TRUE \\
        --cleanUp FALSE \\
        --gatkDir /usr/GenomeAnalysisTK.jar \\
        --novoDir /usr/bin/ \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lohhla: "1.1.7"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.DNA.HLAlossPrediction_CI.xls
    touch ${prefix}_result.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lohhla: "stub"
    END_VERSIONS
    """
}
