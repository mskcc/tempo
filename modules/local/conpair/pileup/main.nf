// CONPAIR Pileup — matches original Tempo QcPileup
// Genome-aware marker bed selection + GATK pileup
process CONPAIR_PILEUP {
    tag "$meta.id"
    label 'process_low'
    container "docker.io/cmopipeline/conpair:v0.3.3"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.pileup"), emit: pileup
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = params.genome ?: 'GRCh37'
    def conpairPath = "/usr/bin/conpair"
    def gatkPath = "/usr/bin/GenomeAnalysisTK.jar"
    def markersBed = genome == 'GRCh38' ?
        "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed" :
        "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed"
    """
    ${conpairPath}/scripts/run_gatk_pileup_for_sample.py \\
        --gatk=${gatkPath} \\
        --bam=${bam} \\
        --markers=${markersBed} \\
        --reference=${fasta} \\
        --outfile=${prefix}.pileup

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "0.3.3"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pileup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        conpair: "stub"
    END_VERSIONS
    """
}
