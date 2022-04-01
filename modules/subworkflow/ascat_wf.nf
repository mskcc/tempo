include { runAscatAlleleCount }		from '../process/Ascat/runAscatAlleleCount' 
include { runAscat }			from '../process/Ascat/runAscat' 

workflow ascat_wf {
	take:
	paired_bams

	main:
	referenceMap = params.referenceMap

	ascatAlleleCountLimit = ["GRCh37","smallGRCh37","GRCh38"].contains(params.genome) ? 48 : 1 
	ascatAlleleCountSegments = Channel.from(1..ascatAlleleCountLimit)
	runAscatAlleleCount(
		ascatAlleleCountSegments,
		ascatAlleleCountLimit,
		paired_bams,
		Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.snpGcCorrections])
	)

	runAscat(
		runAscatAlleleCount.out
		.groupTuple(by:[0,1,2], size: ascatAlleleCountLimit)
		.map{idTumor, idNormal, target, ascatTar ->
			[idTumor, idNormal, target, ascatTar.flatten() ]
		}.combine(paired_bams, by:[0,1,2]),
		Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex, referenceMap.snpGcCorrections])
	)

	emit:
	ascatCNV = runAscat.out.caveman
	ascatSS  = runAscat.out.samplestatistics

}


