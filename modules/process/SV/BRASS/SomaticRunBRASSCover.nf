process runBRASSCover {
  tag "${idTumor}__${idNormal}@${brassCoverIndex}"
  label 'BRASS'

  input:
    each brassCoverIndex
    val(brassCoverLimit)
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(basTumor), path(bamNormal), path(baiNormal), path(basNormal)
    path(genomeFile)
    path(genomeIndex)
	  path("brassRefDir")
	  path(vagrentRefDir)

  output:
  tuple val(idTumor), val(idNormal), val(target), path("brass/tmpBrass/cover/*.*"), path("brass/tmpBrass/progress/*.*")

  script:
  if (params.genome in ["GRCh37","smallGRCh37"]){
    species = "HUMAN"
    assembly = 37   
  }
  else if (params.genome in ["GRCh38"]) {
    species = "HUMAN"
    assembly = 38
  } else { // not sure if run will complete with these params. 
    species = params.genome 
    assembly = params.genome
  }
  if (brassCoverLimit == 1 ) { 
    indexParam = ""
  } else {
    indexParam = "-i ${brassCoverIndex} -l ${brassCoverLimit}"
  }
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass
    for i in rho Ploidy GenderChr GenderChrFound ; do echo \$i ;done > samplestatistics.txt
    brass.pl -j 4 -k 4 -c ${task.cpus} \\
        -d brassRefDir/HiDepth.bed.gz \\
        -f brassRefDir/brass_np.groups.gz \\
        -g ${genomeFile} \\
        -s "${species}" -as "${assembly}" -pr "WGS" \\
        -g_cache ${vagrentRefDir}/vagrent.cache.gz \\
        -vi brassRefDir/viral.genomic.fa.2bit \\
        -mi brassRefDir/all_ncbi_bacteria \\
        -b brassRefDir/500bp_windows.gc.bed.gz \\
        -ct brassRefDir/CentTelo.tsv \\
        -cb brassRefDir/cytoband.txt \\
        -t ${bamTumor} \\
        -n ${bamNormal} \\
        -ss samplestatistics.txt \\
        -o brass \\
        -p cover \\
        ${indexParam}
    """
}

