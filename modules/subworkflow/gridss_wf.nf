#!/usr/bin/env nextflow

// nextflow.enable.dsl=2
include { PreprocessRunGRIDSS }   from '../process/SV/GRIDSS/PreprocessRunGRIDSS'
include { AssemblyRunGRIDSS }     from '../process/SV/GRIDSS/AssemblyRunGRIDSS'
include { VariantCallGRIDSS }     from '../process/SV/GRIDSS/VariantCallGRIDSS'
include { FilterGRIPPS }          from '../process/SV/GRIDSS/FilterGRIPPS'


workflow gridss_wf {
    take:

    bamFiles
    

    main:
    referenceMap = params.referenceMap

    inputPairing = bamFiles
		.map{ row -> 
			[row[0],row[1],row[3],row[5]]
		}

    
    bamFiles
		.map{ [it[0],it[3]]
		}.set {bamTumorFiles}


    bamNormalFiles = bamFiles
		.map{ [it[1],it[5]]
		}
    
    bamFilesUnpaired = bamFiles
		.map{ row -> [row[0],row[3],row[0],row[1]]}
		.mix(
			bamFiles
				.map{ row -> [row[1],row[5],row[0],row[1]]}
				.unique()
		)

    
    PreprocessRunGRIDSS(
		bamFilesUnpaired,
        referenceMap.genomeFile,
        referenceMap.svCallingExcludeRegions,
        referenceMap.bwaIndex,
        referenceMap.genomeDict,
        referenceMap.genomeIndex
	)

    // PreprocessRunGRIDSS.out[0].view()

	bamPairing = inputPairing
		.combine(PreprocessRunGRIDSS.out[0])
		.branch{ idTumor, idNormal,bamTumor,bamNormal, idSample, basename, accessories ->
          tumor: idSample == idTumor
          normal: idSample == idNormal
        }

    tumorNormalPair= bamPairing.tumor
				.combine(bamPairing.normal, by:[0,1,2,3]).unique()

    tumorNormalPair.unique()

    Channel.from(0..3).set{ iterChannel }
    AssemblyRunGRIDSS(
		iterChannel,
        tumorNormalPair,
        referenceMap.genomeFile,
        referenceMap.svCallingExcludeRegions,
        referenceMap.bwaIndex,
        referenceMap.genomeDict,
        referenceMap.genomeIndex
	)

    groupedAssembly = AssemblyRunGRIDSS.out[0].groupTuple(size:4,by:[0,4]).map{it -> [it[1][0],it[2][0],it[3][0],it[4][0],it[5][0],it[6][0],it[7][0],it[7][1],it[7][2],it[7][3]]}  //, val(idTumor), path('mapT'), val(idNormal), path('mapN'),path('assembly1_'),path('assembly2_'),path('assembly3_'),path('assembly4_')


    varCallT = bamTumorFiles.combine(groupedAssembly, by: 0).map{
        idTumor,bamT,baseTumor, mapT, idNormal,baseNormal,mapN,assembly1_,assembly2_,assembly3_,assembly4_ -> [idNormal,bamT,idTumor,baseTumor, mapT, baseNormal,mapN,assembly1_,assembly2_,assembly3_,assembly4_ ]
    }

    varCallIn = varCallT.combine(bamNormalFiles, by: 0).map{
        idNormal,bamT,idTumor,baseTumor, mapT, baseNormal,mapN,assembly1_,assembly2_,assembly3_,assembly4_,bamN -> [bamN,bamT,idTumor,baseTumor, mapT, idNormal,baseNormal,mapN,assembly1_,assembly2_,assembly3_,assembly4_ ]
    }

    VariantCallGRIDSS(
		varCallIn,
        referenceMap.genomeFile,
        referenceMap.svCallingExcludeRegions,
        referenceMap.bwaIndex,
        referenceMap.genomeDict,
        referenceMap.genomeIndex
    )



    // FilterGRIPPS(
    //     VariantCallGRIDSS.out[0],
    //     referenceMap.genomeFile,
    //     referenceMap.svCallingExcludeRegions,
    //     referenceMap.bwaIndex,
    //     referenceMap.genomeDict,
    //     referenceMap.genomeIndex,
    // )



	emit:
	// gridssOutput= FilterGRIPPS.out[1]
    VariantCallGRIDSS.out[0]
}

// workflow {
    
//     bamFiles = Channel.of(["DO6436_tumor", "DO6436_normal", "wgs", "/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_tumor.bam", "/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_tumor.bam.bai", "/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_normal.bam", "/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_normal.bam.bai"] )
//     // bamFiles = Channel.of(['idTumor':"DO6436_tumor", 'idNormal':"DO6436_normal", 'target':"wgs", 'bamTumor':"/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_tumor.bam", 'baiTumor':"/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_tumor.bam.bai", 'bamNormal':"/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_normal.bam", 'baiNormal':"/juno/work/ccs/orgeraj/tempoDev/tstbams/DO6436_normal.bam.bai"] )

//     // tst.map { item -> // re-order the elements
//     //             def idTumor = item[0]
//     //             def idNormal = item[1]
//     //             def target = item[2]
//     //             def bamTumor = item[3]
//     //             def baiTumor = item[4]
//     //             def bamNormal = item[5]
//     //             def baiNormal = item[6]

//     //             return [ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal ]
//     //         }
//     //         .set{ bamFiles }
//     idTumor= "DO6436_tumor"
//     idNormal = "DO6436_normal"
//     gridss_wf( bamFiles )
// }