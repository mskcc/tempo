process CONPAIR_ALL {
    tag "${meta.tumor_id}__${meta.normal_id}"
    label 'process_medium'
    container 'cmopipeline/conpair:v0.3.3'

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
    """
    touch .Rprofile

    echo "${meta.normal_id}\t${meta.tumor_id}" > pairing.txt

    verify_concordances.py \\
        --tumor_pileup=${pileup_tumor} \\
        --normal_pileup=${pileup_normal} \\
        --pairing=pairing.txt \\
        --normal_homozygous_markers_only \\
        --outpre=${prefix}

    estimate_tumor_normal_contaminations.py \\
        --tumor_pileup=${pileup_tumor} \\
        --normal_pileup=${pileup_normal} \\
        --pairing=pairing.txt \\
        --outpre=${prefix}

    mv ${prefix}_concordance.txt ${prefix}.concordance.txt
    mv ${prefix}_contamination.txt ${prefix}.contamination.txt
    """
}
