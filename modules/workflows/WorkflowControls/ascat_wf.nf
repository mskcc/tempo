include { runAscatAlleleCount }		from '../HRDetect/runAscatAlleleCount' 
include { runAscat }				from '../HRDetect/runAscat' 

workflow ascat_wf {
	take:
	paired_bams

	main:
	referenceMap = params.referenceMap

	ascatAlleleCountSegments = Channel.from(1..params.ascat.ascatAlleleCountLimit)
	runAscatAlleleCount(
		ascatAlleleCountSegments,
		params.ascat.ascatAlleleCountLimit,
		paired_bams,
		Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.snpGcCorrections])
	)

	runAscat(
		runAscatAlleleCount.out
		.groupTuple(by:[0,1,2], size: params.ascat.ascatAlleleCountLimit)
		.map{idTumor, idNormal, target, ascatTar ->
			[idTumor, idNormal, target, ascatTar.flatten() ]
		}.combine(paired_bams, by:[0,1]),
		Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.snpGcCorrections])
	)

	emit:
	ascatCNV = runAscat.out

}


