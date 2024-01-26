process VariantCallGRIDSS {
    tag "${idTumor}_${idNormal}_VariantCall"
    publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/GRIDDS", mode: params.publishDirMode, pattern: "*.{vcf.gz,vcf.gz.tbi}"
    
    input:
    tuple path(normalBam), path(tumorBam),val(idTumor),val(baseTumor), path('mapT'), val(idNormal),val(baseNormal), path('mapN'),path('assembly1_'),path('assembly2_'),path('assembly3_'),path('assembly4_')
    path(genomeFile)
    path(excludeRegions)
    path(bwaIndex)
    path(genomeDict)
    path(genomeIndex)

    output:
    tuple val(idTumor), val(idNormal), path("${idTumor}.gripss.filtered*"), emit: GRIDSS4Combine
    path("*.vcf*"), emit: allVcfs
    path('*grids*log') , emit: logs

    
    script:

    """

        grep -v '^chr' ${excludeRegions} > exclude_.bed

        baseTumorBam=`basename $tumorBam`
        baseNormalBam=`basename $normalBam`
        
        mkdir assembly.bam.gridss.working

        mkdir ${baseTumor}.gridss.working

        mkdir ${normalBam}.gridss.working
        
        echo $baseTumor

        for FILE in \$(find . -type l); do
            tgt=`readlink "\${FILE}"`
            filesep=`basename \$tgt `
            echo \$FILE
            echo \$tgt
            echo \$filesep
            if [[ "\$filesep" == assembly* ]]; then
                mv \$FILE ./assembly.bam.gridss.working/\${filesep}

            elif [[ \$filesep == ${baseTumor}* ]]; then
                echo \$filesep
                echo "equals base tumor * "

                if [[ \$filesep == \${baseTumorBam} ]]; then
                    echo ${tumorBam}
                else

                    mv \$FILE ./${baseTumor}.gridss.working/\${filesep}

                fi
            

            elif [[ \$filesep == ${baseNormal}* ]]; then
                if [[ \$filesep == \${baseNormalBam} ]]; then
                    echo ${normalBam}
                else

                    mv \$FILE ./${baseNormal}.gridss.working/\${filesep}
                fi
            fi

        done;



        /opt/gridss/gridss \
        --jvmheap 112g --otherjvmheap 16g \
        -r ${genomeFile} \
        -j /opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
         -t 8 \
         -s assemble,call \
         -b  exclude_.bed \
        -a assembly.bam \
         -o ${idTumor}.GRIDSS.vcf \
        ${normalBam} ${tumorBam} 

        java -Xmx230g -jar /opt/gripss_v2.3.4.jar \
         -sample ${idTumor} \
         -reference ${idNormal} \
         -ref_genome_version 37 \
         -ref_genome ${genomeFile} \
         -known_hotspot_file /juno/work/ccs/orgeraj/SV/sv/known_fusions.37.bedpe \
         -repeat_mask_file /juno/work/ccs/orgeraj/SV/sv/repeat_mask_data.37.fa.gz \
         -vcf ${idTumor}.GRIDSS.vcf \
         -output_dir ./ 

    """

}
