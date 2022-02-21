include { HRDetect }				from '../HRDetect/HRDetect'
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
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.snpGcCorrections
	)

	runAscat(
		runAscatAlleleCount.out
		.groupTuple(by:[0,1,2], size: ascatAlleleCountLimit)
		.map{idTumor, idNormal, target, ascatTar ->
			[idTumor, idNormal, target, ascatTar.flatten() ]
		}.combine(paired_bams, by:[0,1]),
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.snpGcCorrections
	)

	emit:
	ascatCNV = runAscat.out.caveman
	ascatSS = runAscat.out.samplestatistics

}


workflow hrdetect_wf {
	take:
	CNVoutput
	MAFoutput 
	SVoutput
	
	main:
	HRDetectVariantsIn = MAFoutput
		.combine(CNVoutput,by:[0,1,2])
		.combine(SVoutput, by:[0,1,2])

	HRDetect(
		HRDetectVariantsIn,
		workflow.projectDir + "/containers/hrdetect/HRDetect.R"
	)
	
}
