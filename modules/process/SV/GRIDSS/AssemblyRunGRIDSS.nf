process AssemblyRunGRIDSS {
    tag "${idTumor}_${idNormal}"

    input:
    each iterChannel
    tuple val(idTumor), val(idNormal), path(tumorBam),path(normalBam), val(idTumordup1),val(baseTumor), path('liT'), val(idNormaldup1), val(baseNormal),  path('liN')
    path(genomeFile)
    path(excludeRegions)
    path(bwaIndex)
    path(genomeDict)
    path(genomeIndex)

    output:
    tuple val(idTumor), val(idTumordup1),val(baseTumor), path("${baseTumor}.gridss.working/*.bam*"), val(idNormal), val(baseNormal), path("${baseNormal}.gridss.working/*.bam*"),path("assembly.bam.gridss.working/*")
    path "exclude_${iterChannel}.bed"
    path '*gridss*'
    
    script:

    gridds_jar="/opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar"

    """

        echo ${idTumor}
        echo ${tumorBam}

        grep -v '^chr' ${excludeRegions} > exclude_${iterChannel}.bed

        mkdir ${baseNormal}.gridss.working
        mkdir ${baseTumor}.gridss.working

        baseTumorBam=`basename $tumorBam`
        baseNormalBam=`basename $normalBam`

        echo \$baseTumorBam
        echo \$baseNormalBam
        for FILE in \$(find . -type l); do

            tgt=`readlink "\${FILE}"`
            filesep=`basename \$tgt `
            echo \$FILE
            echo \$tgt
            echo \$filesep

            fileDotSep=`echo \$filesep | tr "." "\n"`


            if [[ \$filesep == ${baseTumor}* ]]; then
                if [[ \$filesep == \${baseTumorBam} ]]; then
                    echo ${tumorBam}
                else

                    mv \$FILE ./${baseTumor}.gridss.working/\${filesep}

                fi
            fi

            if [[ \$filesep == ${baseNormal}* ]]; then
                if [[ \$filesep == \${baseNormalBam} ]]; then
                    echo ${normalBam}
                else

                    mv \$FILE ./${baseNormal}.gridss.working/\${filesep}
                fi
            fi

        done

        
        ls
        ls ./${baseNormal}.gridss.working/


        /opt/gridss/gridss \
        --jvmheap ${task.memory.toGiga() - 1}g \
        --otherjvmheap ${task.memory.toGiga() - 1}g \
        -r ${genomeFile} \
        -j ${gridds_jar} \
        -t 8 \
        -s assemble \
        -a assembly.bam \
        -b exclude_${iterChannel}.bed \
        --jobnodes 4 \
        --jobindex ${iterChannel} \
        --workingdir ./ \
        ${normalBam} ${tumorBam}  || cat *gridss*.log 


    """

}
