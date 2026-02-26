process SVCLONE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/svclone:0.0.1' :
        'cmopipeline/svclone:0.0.1' }"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(bedpe), path(maf), path(cnv), path(ploidy)

    output:
    tuple val(meta), path("${prefix}"),                          emit: output_dir
    tuple val(meta), path("svclone/svs/*cluster_certainty.txt"), path("svclone/snvs/*cluster_certainty.txt"), emit: cluster_certainty
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /opt/svclone/svclone_wrapper.py \\
        --cfg_template /config/svclone_config.ini \\
        --bedpe ${bedpe} \\
        --maf ${maf} \\
        --purity_ploidy ${ploidy} \\
        --out_dir svclone_in \\
        --sampleid ${prefix} \\
        --bam ${tumor_bam} \\
        --cnv ${cnv}

    mkdir -p svclone/svs svclone/snvs
    cp ${prefix}/ccube_out/post_assign/*.RData ${prefix}/ccube_out/post_assign/*.pdf svclone/svs/ || true
    cp ${prefix}/ccube_out/post_assign/snvs/*.RData ${prefix}/ccube_out/post_assign/snvs/*.pdf svclone/snvs/ || true
    for i in ${prefix}/ccube_out/post_assign/*.txt ; do
        sed "s/^/${prefix}\\t/g" \$i | sed "0,/^${prefix}\\t/s//sampleid\\t/" > svclone/svs/\$(basename \$i)
    done
    for i in ${prefix}/ccube_out/post_assign/snvs/*.txt ; do
        sed "s/^/${prefix}\\t/g" \$i | sed "0,/^${prefix}\\t/s//sampleid\\t/" > svclone/snvs/\$(basename \$i)
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svclone: 0.0.1
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix} svclone/svs svclone/snvs
    touch svclone/svs/sv_stub_cluster_certainty.txt svclone/snvs/snv_stub_cluster_certainty.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svclone: 0.0.1
    END_VERSIONS
    """
}
