// CONPAIR All — matches original Tempo QcConpairAll
// Genome-aware marker txt selection + pairing file
process CONPAIR_ALL {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'docker.io/cmopipeline/conpair:v0.3.3'

    input:
    tuple val(meta), path(pileup_tumor), path(pileup_normal)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("*.concordance.txt"), path("*.contamination.txt"), emit: conpair_output

    stub:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    """
    touch ${prefix}.concordance.txt ${prefix}.contamination.txt
    """

    script:
    def prefix = "${meta.tumor_id}__${meta.normal_id}"
    def genome = params.genome ?: 'GRCh37'
    def conpairPath = "/usr/bin/conpair"
    def markersTxt = genome == 'GRCh38' ?
        "${conpairPath}/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt" :
        "${conpairPath}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt"
    """
    touch .Rprofile

    echo "${meta.normal_id}\\t${meta.tumor_id}" > pairing.txt

    ${conpairPath}/scripts/verify_concordances.py \\
        --tumor_pileup=${pileup_tumor} \\
        --normal_pileup=${pileup_normal} \\
        --markers=${markersTxt} \\
        --pairing=pairing.txt \\
        --normal_homozygous_markers_only \\
        --outpre=${prefix}

    ${conpairPath}/scripts/estimate_tumor_normal_contaminations.py \\
        --tumor_pileup=${pileup_tumor} \\
        --normal_pileup=${pileup_normal} \\
        --markers=${markersTxt} \\
        --pairing=pairing.txt \\
        --outpre=${prefix}

    mv ${prefix}_concordance.txt ${prefix}.concordance.txt
    mv ${prefix}_contamination.txt ${prefix}.contamination.txt
    """
}
