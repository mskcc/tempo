process PreprocessRunGRIDSS {
    tag "${idSample}"
    
    input:
	tuple val(idSample), path(bamPath),val(idTumor),val(idNormal)
	path(genomeFile)
	path(excludeRegions)
    path(bwaIndex)
    path(genomeDict)
    path(genomeIndex)

	output:
    tuple val(idSample), env(baseBam), path("*.gridss.working/*.bam*")
    path '*gridss*'


	script:
	gripss_jar="/opt/hmftools/gripss-1.9.jar"
	gRIDSS_jar="/opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar"



    """

        /opt/gridss/gridss \
        -r ${genomeFile} \
        -j ${gRIDSS_jar} \
        -s preprocess \
        -t 4 \
        -b ${excludeRegions} \
        ${bamPath}

        cigarbam="${idSample}.bam.cigar_metrics"
        samtags="${idSample}.bam.computesamtags.changes.tsv"
        blacklistcov="${idSample}.bam.coverage.blacklist.bed"
        idsv_metrics="${idSample}.bam.idsv_metrics"
        insert_metrics="${idSample}.bam.insert_size_metrics"
        mapqmetrics="${idSample}.bam.mapq_metrics"
        SVbam="${idSample}.bam.sv.bam"
        SVbamcsi="${idSample}.bam.sv.bam.csi"
        tagmetrics="${idSample}.bam.tag_metrics"
        baseBam=`basename ${bamPath}`


    """

}

