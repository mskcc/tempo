include { 
	generateBasFile;
	runBRASSInput;
	runBRASSCover;
	runBRASS } from '../SV/BRASS/SomaticRunBRASS.nf'

workflow brass_wf
{
	take:
	bamFiles
	bamFilesPaired
	inputPairing
	sampleStatistics // from ascat 

	main:
	generateBasFile(bamFiles)

	basPairing = inputPairing
		.combine(generateBasFile.out)
		.branch{ idTumor, idNormal, idSample, target, basFile ->
          tumor: idSample == idTumor
          normal: idSample == idNormal
        }

	brassInfiles = bamFilesPaired.combine(basPairing.tumor
		.combine(basPairing.normal, by:[0,1])
		.map{ idTumor,idNormal, idSample1, target1, basFile1, idSample2, target2, c ->
          [idTumor,idNormal,target1,basFile1,basFile2]
        }, by:[0,1,2])

	BRASSInputSegments = Channel.from(1..2)
	runBRASSInput(
		BRASSInputSegments,
		brassInfiles,
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.brassRefDir, 
		referenceMap.vagrentRefDir
		)

	brassCoverLimit = params.genome in ["GRCh37","smallGRCh37","GRCh37"] ? 24 : 1
	BRASSCoverSegments = Channel.from(1..brassCoverLimit)
	runBRASSCover(
		BRASSCoverSegments,
		brassInfiles,
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.brassRefDir, 
		referenceMap.vagrentRefDir
		)

	brassInfilesWithPreprocess = brassInfiles
		.combine(runBRASSInput.out, by:[0,1,2])
		.combine(runBRASSCover.out, by:[0,1,2])
		.combine(sampleStatistics, by:[0,1,2])

	runBRASS(
		brassInfilesWithPreprocess,
		referenceMap.genomeFile, 
		referenceMap.genomeIndex, 
		referenceMap.brassRefDir, 
		referenceMap.vagrentRefDir
	)

	emit:
	brassOutput = runBRASS.out

}