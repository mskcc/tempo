process runBRASS {
    tag "${idTumor}__${idNormal}"
    label 'BRASS'

    publishDir "${params.outDir}/somatic/${outputPrefix}/", mode: params.publishDirMode, pattern: "brass/*.{gz,tbi}"

    input:
    tuple val(idTumor), val(idNormal), val(target), 
      path(bamTumor), path(baiTumor), path(basTumor), 
      path(bamNormal), path(baiNormal), path(basNormal), 
      path(BrassInputTmp), 
      path(BrassInputProgress), 
      path(BrassCoverTmp), 
      path(BrassCoverProgress), 
      path(ascatSampleStatistics) 
    path(genomeFile)
    path(genomeIndex)
    path("brassRefDir")
    path(vagrentRefDir)

    output:
    tuple val(idTumor), val(idNormal), val(target), path("brass/*.{vcf.gz,vcf.gz.tbi}"), emit: BRASSOutput
    tuple val(idTumor), val(idNormal), val(target), path("brass/*.vcf.gz"), path("brass/*.vcf.gz.tbi"), emit: BRASS4Combine


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
        -ss ${ascatSampleStatistics} \\
        -o brass
    """
}
