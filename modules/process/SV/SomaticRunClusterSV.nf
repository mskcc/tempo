process SomaticRunClusterSV {
  tag { outputPrefix }

  publishDir "${params.outDir}/somatic/${outputPrefix}/sv_clustering", mode: params.publishDirMode, pattern: "out/*.tsv"
  
  input:
    tuple val(idTumor), val(idNormal), val(target), path(bedpe)
    
  output:
    tuple val(idTumor), val(idNormal), val(target), path("out/${outputPrefix}.sv_clusters_and_footprints.tsv"), path("${outputPrefix}.sv_distance_pvals")

  when: ["GRCh37","GRCh38","smallGRCh37"].contains(params.genome)
    
  script:
  outputPrefix = [idTumor,idNormal].unique()
  outputPrefix.remove("")
  outputPrefix = outputPrefix.join("__")
  genome_ = params.genome == "GRCh37" || params.genome == 'smallGRCh37' ? "hs37d5" : "hg38"
  """
  mkdir -p tmp
  grep -v "^#" ${bedpe} | cut -f 1-10 > tmp/${outputPrefix}.bedpe
  Rscript /opt/ClusterSV/R/run_cluster_sv.R \\
    -chr /opt/ClusterSV/references/${genome_}.chrom_sizes \\
    -cen_telo /opt/ClusterSV/references/${genome_}_centromere_and_telomere_coords.txt \
    -out ${outputPrefix} \
    -bedpe tmp/${outputPrefix}.bedpe

  mkdir -p out
  grep "^#CHROM" ${bedpe} | cut -f 1-10 | tr -d "\\n" > out/${outputPrefix}.sv_clusters_and_footprints.tsv
  echo -e "\\tcluster_id\\tcluster_total_count\\tfootprint_id_low\\tfootprint_id_high\\tcoord_footprint_id_low\\tcoord_footprint_id_high\\tpval" >> out/${outputPrefix}.sv_clusters_and_footprints.tsv 
  cat ${outputPrefix}.sv_clusters_and_footprints.tsv >> out/${outputPrefix}.sv_clusters_and_footprints.tsv
  """
}
