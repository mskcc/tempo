process SomaticAnnotateSVBedpe {
  tag "${idTumor}__${idNormal}"

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.unfiltered.bedpe"
  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs", mode: params.publishDirMode, pattern: "*.final.bedpe"


  input:
    tuple val(idTumor), val(idNormal), val(target), path(bedpein) 
    path(repeatMasker)
    path(mapabilityBlacklist)
    path(svBlacklistBed)
    path(svBlacklistBedpe)
    path(svBlacklistFoldbackBedpe)
    path(svBlacklistTEBedpe)
    path(spliceSites)
    path(custom_scripts) 
    val(genome)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.unfiltered.bedpe"), emit: SVAnnotBedpe
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.final.bedpe"), emit: SVAnnotBedpePass
    tuple val(idTumor), val(idNormal), path("${outputPrefix}.final.bedpe"), emit: SVAnnotBedpe4Aggregate

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  genome_ = ["GRCh37","smallGRCh37"].contains(genome) ? "hg19" : genome == "GRCh38" ? "hg38" : "hg18"
  """
  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${mapabilityBlacklist} \\
    --bedpe ${bedpein} \\
    --tag mappability \\
    --output ${outputPrefix}.combined.dac.bedpe \\
    --match-type either 
  
  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${repeatMasker} \\
    --bedpe ${outputPrefix}.combined.dac.bedpe \\
    --tag repeat_masker \\
    --output ${outputPrefix}.combined.dac.rm.bedpe \\
    --match-type either

  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${svBlacklistBed} \\
    --bedpe ${outputPrefix}.combined.dac.rm.bedpe \\
    --tag pcawg_blacklist_bed \\
    --output ${outputPrefix}.combined.dac.rm.pcawg.1.bedpe \\
    --match-type either

  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${svBlacklistBedpe} \\
    --bedpe ${outputPrefix}.combined.dac.rm.pcawg.1.bedpe \\
    --tag pcawg_blacklist_bedpe \\
    --output ${outputPrefix}.combined.dac.rm.pcawg.2.bedpe \\
    --match-type both \\
    --ignore-strand

  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${svBlacklistFoldbackBedpe} \\
    --bedpe ${outputPrefix}.combined.dac.rm.pcawg.2.bedpe \\
    --tag pcawg_blacklist_fb_bedpe \\
    --output ${outputPrefix}.combined.dac.rm.pcawg.3.bedpe \\
    --match-type both
  
  python ${custom_scripts}/filter_regions_bedpe.py \\
    --blacklist-regions ${svBlacklistTEBedpe} \\
    --bedpe ${outputPrefix}.combined.dac.rm.pcawg.3.bedpe \\
    --tag pcawg_blacklist_te_bedpe \\
    --output ${outputPrefix}.combined.dac.rm.pcawg.4.bedpe \\
    --match-type either

  python ${custom_scripts}/detect_cdna.py \\
    --exon-junct ${spliceSites} \\
    --bedpe ${outputPrefix}.combined.dac.rm.pcawg.4.bedpe \\
    --out-bedpe ${outputPrefix}.combined.dac.rm.pcawg.cdna.bedpe \\
    --out ${outputPrefix}.contamination.tsv

  python ${custom_scripts}/run_iannotatesv.py \\
    --bedpe ${outputPrefix}.combined.dac.rm.pcawg.cdna.bedpe \\
    --genome ${genome_} \\
    --threads ${task.cpus * 2}

  cp ${outputPrefix}.combined.dac.rm.pcawg.cdna.iannotate.bedpe \\
    ${outputPrefix}.unfiltered.bedpe

  awk -F"\\t" '\$1 ~ /#/ || \$12 == "PASS"' \\
    ${outputPrefix}.unfiltered.bedpe > \\
    ${outputPrefix}.final.bedpe
  
  """

}
