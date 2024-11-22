process SomaticRunSVCircos {
	tag "${idTumor}__${idNormal}"

	publishDir "${params.outDir}/somatic/${outputPrefix}/combined_svs/", mode: params.publishDirMode

	input:
	  tuple val(idTumor), val(idNormal), val(target), path(bedpe), path(cnv)
	  path(biocircos_script)
	  path(biocircos_Rmd)

	output:
	  path("${outputPrefix}.circos.html")

	when: ["GRCh37","smallGRCh37","GRCh38"].contains(params.genome)

	script:
	outputPrefix = "${idTumor}__${idNormal}"
	genome_version = params.genome == 'GRCh38' ? "hg38" : "hg19"
	"""
	Rscript ${biocircos_script} \\
	  -b ${bedpe} \\
	  -c ${cnv} \\
	  -s ${outputPrefix} \\
	  -g ${genome_version}

	"""
}
