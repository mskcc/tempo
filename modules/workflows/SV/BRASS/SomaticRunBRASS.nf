process generateBasFile {
  tag { idSample }

  input:
  set val(idSample), val(target), file(bam), file(bai) 
  path(genomeFile)
  path(genomeIndex)

  output: 
  tuple val(idSample), val(target), path("*.bas")

  script:
  """
  bam_stats \\
    -i ${bam} \\
    -o ${bam}.bas \\
    -r ${genomeIndex} \\
    -@ ${ task.cpus > 1 ? task.cpus - 1 : task.cpus }
  """
}

process runBRASSInput {
    tag { idTumor + "__" + idNormal + "@" + inputIndex }
    label 'BRASS'   

    input:
    each inputIndex
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(basTumor), path(bamNormal), path(baiNormal), path(basNormal)
    path(genomeFile)
    path(genomeIndex)
	  path(brassRefDir)
	  path(vagrentRefDir)

    output:
    tuple val(idTumor), val(idNormal), val(target), path("brass/tmpBrass/*.*"), path("brass/tmpBrass/progress/*.*")

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
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass
    for i in rho Ploidy GenderChr GenderChrFound ; do echo \$i ;done > samplestatistics.txt
    brass.pl -j 4 -k 4 -c ${task.cpus} \\
        -d ${brassRefDir}/HiDepth.bed.gz \\
        -f ${brassRefDir}/brass_np.groups.gz \\
        -g ${genomeFile} \\
        -s "${species}" -as "${assembly}" -pr "WGS" \\
        -g_cache ${vagrentRefDir}/vagrent.cache.gz \\
        -vi ${brassRefDir}/viral.genomic.fa.2bit \\
        -mi ${brassRefDir}/all_ncbi_bacteria \\
        -b ${brassRefDir}/500bp_windows.gc.bed.gz \\
        -ct ${brassRefDir}/CentTelo.tsv \\
        -cb ${brassRefDir}/cytoband.txt \\
        -t ${bamTumor} \\
        -n ${bamNormal} \\
        -ss samplestatistics.txt \\
        -o brass \\
        -p input \\
        -i ${inputIndex} -l 2
    """  
}

process runBRASSCover {
  tag { idTumor + "__" + idNormal + "@" + coverIndex }
  label 'BRASS'

  input:
    each coverIndex
    val(coverLimit)
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(basTumor), path(bamNormal), path(baiNormal), path(basNormal)
    path(genomeFile)
    path(genomeIndex)
	  path(brassRefDir)
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
    indexParam = "-i ${coverIndex} -l ${brassCoverLimit}"
  }
    """
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass
    for i in rho Ploidy GenderChr GenderChrFound ; do echo \$i ;done > samplestatistics.txt
    brass.pl -j 4 -k 4 -c ${task.cpus} \\
        -d ${brassRefDir}/HiDepth.bed.gz \\
        -f ${brassRefDir}/brass_np.groups.gz \\
        -g ${genomeFile} \\
        -s "${species}" -as "${assembly}" -pr "WGS" \\
        -g_cache ${vagrentRefDir}/vagrent.cache.gz \\
        -vi ${brassRefDir}/viral.genomic.fa.2bit \\
        -mi ${brassRefDir}/all_ncbi_bacteria \\
        -b ${brassRefDir}/500bp_windows.gc.bed.gz \\
        -ct ${brassRefDir}/CentTelo.tsv \\
        -cb ${brassRefDir}/cytoband.txt \\
        -t ${bamTumor} \\
        -n ${bamNormal} \\
        -ss samplestatistics.txt \\
        -o brass \\
        -p cover \\
        ${indexParam}
    """
}

process runBRASS {
    tag { idTumor + "__" + idNormal }
    label 'BRASS'

    publishDir "${outDir}/somatic/${outputPrefix}/", mode: params.publishDirMode, pattern: "brass/*.{gz,tbi}"

    input:
    tuple val(idTumor), val(idNormal), val(target), path(bamTumor), path(baiTumor), path(basTumor), path(bamNormal), path(baiNormal), path(basNormal), file(BrassInputTmp), file(BrassInputProgress), file(BrassCoverTmp), file(BrassCoverProgress), file(ascatSampleStatistics) from bamsForBRASSSV
    path(genomeFile)
    path(genomeIndex)
	  path(brassRefDir)
  	path(vagrentRefDir)

    output:
    set idTumor, idNormal, target, file("brass/*.{vcf.gz,vcf.gz.tbi}")

    script:
    outputPrefix = "${idTumor}__${idNormal}"
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
    """
    BrassInputResults=( ${BrassInputTmp.join(" ")} )
    BrassInputProgress=( ${BrassInputProgress.join(" ")} )
    BrassCoverCover=( ${BrassCoverTmp.join(" ")} )
    BrassCoverProgress=( ${BrassCoverProgress.join(" ")} )
  
    export TMPDIR=\$(pwd)/tmp ; mkdir -p \$TMPDIR brass/tmpBrass/progress brass/tmpBrass/cover 
    for i in "\${BrassCoverCover[@]}" ; do 
        mv \$i brass/tmpBrass/cover
    done
    for i in "\${BrassCoverProgress[@]}" ; do 
        mv \$i brass/tmpBrass/progress
    done
    for i in "\${BrassInputResults[@]}" ; do 
        mv \$i brass/tmpBrass
    done
    for i in "\${BrassInputProgress[@]}" ; do 
        mv \$i brass/tmpBrass/progress
    done
    brass.pl -j 4 -k 4 -c ${task.cpus} \\
        -d ${brassRefDir}/HiDepth.bed.gz \\
        -f ${brassRefDir}/brass_np.groups.gz \\
        -g ${genomeFile} \\
        -s "${species}" -as "${assembly}" -pr "WGS" \\
        -g_cache ${vagrentRefDir}/vagrent.cache.gz \\
        -vi ${brassRefDir}/viral.genomic.fa.2bit \\
        -mi ${brassRefDir}/all_ncbi_bacteria \\
        -b ${brassRefDir}/500bp_windows.gc.bed.gz \\
        -ct ${brassRefDir}/CentTelo.tsv \\
        -cb ${brassRefDir}/cytoband.txt \\
        -t ${bamTumor} \\
        -n ${bamNormal} \\
        -ss ${ascatSampleStatistics} \\
        -o brass
    """

}
