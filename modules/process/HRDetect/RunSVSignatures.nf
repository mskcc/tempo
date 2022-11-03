process RunSVSignatures {
  tag {idTumor + "__" + idNormal}

  publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/", mode: params.publishDirMode

  input:
    tuple val(idTumor), val(idNormal), val(target), path(bedpe)
    path(sv_signature_script) 

  output:
    tuple val("placeholder"), val(idTumor), val(idNormal),
      path("${outputPrefix}_catalogues.pdf"), path("${outputPrefix}_exposures.tsv")

  script:
  outputPrefix = "${idTumor}__${idNormal}"
  genome_version = params.genome == 'GRCh38' ? "hg38" : "hg19"
  """
  Rscript ${sv_signature_script} \\
    -i ${bedpe} \\
    -g ${genome_version} \\
    -n ${task.cpus} \\
    -s ${outputPrefix}
  """

}
