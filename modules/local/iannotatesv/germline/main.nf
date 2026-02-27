process IANNOTATESV_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/cmopipeline/iannotatesv:0.0.2' :
        'docker.io/cmopipeline/iannotatesv:0.0.2' }"

    input:
    tuple val(meta), path(bedpe)
    path(repeat_masker)
    path(mapability_blacklist)
    path(sv_blacklist_bed)
    path(sv_blacklist_bedpe)
    path(sv_blacklist_foldback_bedpe)
    path(sv_blacklist_te_bedpe)
    path(splice_sites)
    val(genome)

    output:
    tuple val(meta), path("${prefix}.unfiltered.bedpe"), emit: bedpe_unfiltered
    tuple val(meta), path("${prefix}.final.bedpe"),      emit: bedpe_pass
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def genome_ = genome == 'GRCh37' ? 'hg19' : (genome == 'GRCh38' ? 'hg38' : genome)
    """
    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${mapability_blacklist} \\
        --bedpe ${bedpe} --tag mappability \\
        --output ${prefix}.dac.bedpe --match-type either

    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${repeat_masker} \\
        --bedpe ${prefix}.dac.bedpe --tag repeat_masker \\
        --output ${prefix}.dac.rm.bedpe --match-type either

    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${sv_blacklist_bed} \\
        --bedpe ${prefix}.dac.rm.bedpe --tag pcawg_blacklist_bed \\
        --output ${prefix}.dac.rm.pcawg.1.bedpe --match-type either

    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${sv_blacklist_bedpe} \\
        --bedpe ${prefix}.dac.rm.pcawg.1.bedpe --tag pcawg_blacklist_bedpe \\
        --output ${prefix}.dac.rm.pcawg.2.bedpe --match-type both --ignore-strand

    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${sv_blacklist_foldback_bedpe} \\
        --bedpe ${prefix}.dac.rm.pcawg.2.bedpe --tag pcawg_blacklist_fb_bedpe \\
        --output ${prefix}.dac.rm.pcawg.3.bedpe --match-type both

    python /opt/iannotatesv/filter_regions_bedpe.py \\
        --blacklist-regions ${sv_blacklist_te_bedpe} \\
        --bedpe ${prefix}.dac.rm.pcawg.3.bedpe --tag pcawg_blacklist_te_bedpe \\
        --output ${prefix}.dac.rm.pcawg.4.bedpe --match-type either

    python /opt/iannotatesv/detect_cdna.py \\
        --exon-junct ${splice_sites} \\
        --bedpe ${prefix}.dac.rm.pcawg.4.bedpe \\
        --out-bedpe ${prefix}.dac.rm.pcawg.cdna.bedpe \\
        --out ${prefix}.contamination.tsv

    python /opt/iannotatesv/run_iannotatesv.py \\
        --bedpe ${prefix}.dac.rm.pcawg.cdna.bedpe \\
        --genome ${genome_} \\
        --threads ${task.cpus * 2}

    cp ${prefix}.dac.rm.pcawg.cdna.iannotate.bedpe ${prefix}.unfiltered.bedpe
    awk -F"\\t" '\$1 ~ /#/ || \$12 == "PASS"' ${prefix}.unfiltered.bedpe > ${prefix}.final.bedpe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iannotatesv: 0.0.2
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.unfiltered.bedpe ${prefix}.final.bedpe
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iannotatesv: 0.0.2
    END_VERSIONS
    """
}
