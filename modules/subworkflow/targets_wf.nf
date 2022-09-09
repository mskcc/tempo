include { CreateBaitsetFiles }   from '../process/Targets/CreateBaitsetFiles'

workflow targets_wf
{
  main:
    referenceMap = params.referenceMap
    targetsMap   = params.targetsMap

    CreateBaitsetFiles(
    	Channel.from(targetsMap.keySet())
		.map{ targetId ->
			[targetId, targetsMap."${targetId}".targetsBed, targetsMap."${targetId}".baitsBed]
		},
	referenceMap.genomeFile,
	referenceMap.genomeIndex,
	referenceMap.genomeDict,
	referenceMap.codingRegions
    )

  emit:
    baitsetInterval = CreateBaitsetFiles.out.baitsetInterval
    codingBaitsetBed = CreateBaitsetFiles.out.codingBaitsetBed
    baitsetPlus5 = CreateBaitsetFiles.out.baitsetPlus5
    baitsetPlus5_unzipped = CreateBaitsetFiles.out.baitsetPlus5_unzipped
}
