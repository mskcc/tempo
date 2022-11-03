process SomaticRunClusterSV {
  tag { outputPrefix }

  publishDir "${params.outDir}/somatic/${outputPrefix}/sv_clustering", mode: params.publishDirMode, pattern: "*.tsv"
  
  input:
    tuple val(idTumor), val(idNormal), val(target), path(bedpe)
    
  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.sv_clusters_and_footprints.tsv"), path("${outputPrefix}.sv_distance_pvals")

  when: ["GRCh37","GRCh38","smallGRCh37"].contains(params.genome)
    
  script:
  outputPrefix = [idTumor,idNormal].unique()
  outputPrefix.remove("")
  outputPrefix = outputPrefix.join("__")
  genome_ = params.genome == "GRCh37" || params.genome == 'smallGRCh37' ? "hs37d5" : "hg38"
  """
  mkdir -p tmp
  grep -v "^#" ${bedpe} | cut -f 1-12 > tmp/${outputPrefix}.bedpe
  if [ \$(cat ${outputPrefix}.bedpe | wc -l ) -lt 1 ] ; then
    touch tmp/${outputPrefix}.sv_clusters_and_footprints.tsv
  else
    Rscript /opt/ClusterSV/R/run_cluster_sv.R \\
      -chr /opt/ClusterSV/references/${genome_}.chrom_sizes \\
      -cen_telo /opt/ClusterSV/references/${genome_}_centromere_and_telomere_coords.txt \
      -out ${outputPrefix} \
      -bedpe tmp/${outputPrefix}.bedpe

    mv ${outputPrefix}.sv_clusters_and_footprints.tsv tmp/${outputPrefix}.sv_clusters_and_footprints.tsv
  fi

  grep "^#CHROM" ${bedpe} | cut -f 1-12 | tr -d "\\n" > ${outputPrefix}.sv_clusters_and_footprints.tsv
  echo -e "\\tcluster_id\\tcluster_total_count\\tfootprint_id_low\\tfootprint_id_high\\tcoord_footprint_id_low\\tcoord_footprint_id_high\\tpval" \\
    >> ${outputPrefix}.sv_clusters_and_footprints.tsv
  cat tmp/${outputPrefix}.sv_clusters_and_footprints.tsv >> ${outputPrefix}.sv_clusters_and_footprints.tsv
  """
}
