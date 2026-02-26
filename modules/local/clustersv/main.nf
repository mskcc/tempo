process CLUSTERSV {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://cmopipeline/clustersv:0.0.1' :
        'cmopipeline/clustersv:0.0.1' }"

    input:
    tuple val(meta), path(bedpe)
    val(genome)

    output:
    tuple val(meta), path("${prefix}.sv_clusters_and_footprints.tsv"), path("${prefix}.sv_distance_pvals"), path("${prefix}.clustered.bedpe"), emit: clustered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def genome_ = (genome == 'GRCh37' || genome == 'smallGRCh37') ? 'hs37d5' : 'hg38'
    """
    mkdir -p tmp
    grep -v "^#" ${bedpe} | cut -f 1-10 > tmp/${prefix}.bedpe
    if [ \$(cat tmp/${prefix}.bedpe | wc -l) -lt 1 ] ; then
        touch ${prefix}.sv_clusters_and_footprints.tsv
        touch ${prefix}.sv_distance_pvals
    else
        Rscript /opt/ClusterSV/R/run_cluster_sv.R \\
            -chr /opt/ClusterSV/references/${genome_}.chrom_sizes \\
            -cen_telo /opt/ClusterSV/references/${genome_}_centromere_and_telomere_coords.txt \\
            -out ${prefix} \\
            -bedpe tmp/${prefix}.bedpe
    fi

    grep "^##" ${bedpe} > ${prefix}.clustered.bedpe
    grep "^#CHROM" ${bedpe} | tr "\\n" "\\t" >> ${prefix}.clustered.bedpe
    echo -e "cluster_id\\tcluster_total_count\\tfootprint_id_low\\tfootprint_id_high\\tcoord_footprint_id_low\\tcoord_footprint_id_high\\tclustersv_pval" >> ${prefix}.clustered.bedpe
    paste <(grep -v "^#" ${bedpe}) <(cut -f 11- ${prefix}.sv_clusters_and_footprints.tsv) >> ${prefix}.clustered.bedpe || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustersv: 0.0.1
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sv_clusters_and_footprints.tsv ${prefix}.sv_distance_pvals ${prefix}.clustered.bedpe
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clustersv: 0.0.1
    END_VERSIONS
    """
}
