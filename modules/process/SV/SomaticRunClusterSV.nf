process SomaticRunClusterSV {
  tag { outputPrefix }

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.clustered.bedpe"
  
  input:
    tuple val(idTumor), val(idNormal), val(target), path(bedpe)
    
  output:
    tuple val(idTumor), val(idNormal), val(target),
      path("${outputPrefix}.sv_clusters_and_footprints.tsv"),
      path("${outputPrefix}.sv_distance_pvals"),
      path("${clusteredBedpe}"), emit: clusterSvOutput
    tuple val("placeholder"), val(idTumor), val(idNormal),
      path("${clusteredBedpe}"), emit: Bedpe4Aggregate

  when: ["GRCh37","GRCh38","smallGRCh37"].contains(params.genome)
    
  script:
  outputPrefix = [idTumor,idNormal].unique()
  outputPrefix.remove("")
  outputPrefix = outputPrefix.join("__")
  genome_ = params.genome == "GRCh37" || params.genome == 'smallGRCh37' ? "hs37d5" : "hg38"
  clusteredBedpe = bedpe.getBaseName() + ".clustered.bedpe"
  """
  mkdir -p tmp
  grep -v "^#" ${bedpe} | cut -f 1-10 > tmp/${outputPrefix}.bedpe
  if [ \$(cat tmp/${outputPrefix}.bedpe | wc -l ) -lt 1 ] ; then
    touch tmp/${outputPrefix}.sv_clusters_and_footprints.tsv
  else
    Rscript /opt/ClusterSV/R/run_cluster_sv.R \\
      -chr /opt/ClusterSV/references/${genome_}.chrom_sizes \\
      -cen_telo /opt/ClusterSV/references/${genome_}_centromere_and_telomere_coords.txt \
      -out ${outputPrefix} \
      -bedpe tmp/${outputPrefix}.bedpe
  fi

  grep "^##" ${bedpe} > ${clusteredBedpe}
  grep "^#CHROM" ${bedpe} | tr "\\n" "\\t" >> ${clusteredBedpe}
  echo -e "cluster_id\\tcluster_total_count\\tfootprint_id_low\\tfootprint_id_high\\tcoord_footprint_id_low\\tcoord_footprint_id_high\\tclustersv_pval" \\
    >> ${clusteredBedpe}
  paste <(grep -v "^#" ${bedpe} ) \
    <( cut -f 11- ${outputPrefix}.sv_clusters_and_footprints.tsv) >> ${clusteredBedpe}
  """
}
