
workflow PairTumorNormal
{
    take:
        inputBam
        inputPairing

    main:
        if (params.pairing) {
            // Parse input FASTQ mapping
            inputBam.combine(inputPairing)
            .filter { item ->
                def idSample = item[0]
                def target = item[1]
                def sampleBam = item[2]
                def sampleBai = item[3]
                def idTumor = item[4]
                def idNormal = item[5]
                idSample == idTumor
            }.map { item ->
                def idTumor = item[4]
                def idNormal = item[5]
                def tumorBam = item[2]
                def tumorBai = item[3]
                def target = item[1]
                return [ idTumor, idNormal, target, tumorBam, tumorBai ]
            }
            .unique()
            .set{ bamsTumor }

            inputBam.combine(inputPairing)
            .filter { item ->
                def idSample = item[0]
                def target = item[1]
                def sampleBam = item[2]
                def sampleBai = item[3]
                def idTumor = item[4]
                def idNormal = item[5]
                idSample == idNormal
            }.map { item ->
                def idTumor = item[4]
                def idNormal = item[5]
                def normalBam = item[2]
                def normalBai = item[3]
                def target = item[1]
                return [ idTumor, idNormal, target, normalBam, normalBai ]
            }
            .unique()
            .set{ bamsNormal }

            bamsNormal.map { item ->
            def idNormal = item[1]
            def target = item[2]
            def normalBam = item[3]
            def normalBai = item[4]
            return [ idNormal, target, normalBam, normalBai ] }
            .unique()
            .set{ bams }

            bamsTumor.combine(bamsNormal, by: [0,1,2])
            .map { item -> // re-order the elements
                def idTumor = item[0]
                def idNormal = item[1]
                def target = item[2]
                def bamTumor = item[3]
                def baiTumor = item[4]
                def bamNormal = item[5]
                def baiNormal = item[6]

                return [ idTumor, idNormal, target, bamTumor, baiTumor, bamNormal, baiNormal ]
            }
            .set{ bamFiles }
        }
        else
        {
            println "Error: PairTumorNormal - pairing file not provided."
            exit 1
        } 


    emit:
        bamFiles   = bamFiles
        bams       = bams
        bamsNormal = bamsNormal
        bamsTumor  = bamsTumor
}
