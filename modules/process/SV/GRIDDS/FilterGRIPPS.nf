process FilterGRIPPS {
    tag "${idTumor}_${idNormal}"

    input:
    tuple val(idTumor), val(idNormal), path(bamTumor),  path(bamNormal)
    path(genomeFile)
    path(excludeRegions)
    path(bwaIndex)
    path(genomeDict)
    path(genomeIndex)

    output:
    
    tuple val(idTumor), val(idNormal)
    path "${workdir}/${idTumor}_${idNormal}.GRIDSS.grippsfiltered.vcf"
    
    script:

    def gripss_jar="/juno/work/ccs/orgeraj/SV/gripss-1.9.jar"
    def gRIDSS_jar="/opt/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar"
    def workdir="/juno/work/ccs/orgeraj/SV/workdirs/${idTumor}_${idNormal}"

    """

    
        echo $idTumor
        grep -v '^chr' ${excludeRegions} > exclude_.bed


        java -jar ${gripss_jar} \
         -sample ${bamTumor} \
         -reference ${bamNormal} \
         -ref_genome_version 37 \
         -ref_genome ${genomeFile} \
         -b  exclude_.bed \
         -known_hotspot_file /juno/work/ccs/orgeraj/SV/sv/known_fusions.37.bedpe \
         -repeat_mask_file /juno/work/ccs/orgeraj/SV/sv/repeat_mask_data.37.fa.gz \
         -vcf ${workdir}/${idTumor}_${idNormal}.GRIDSS.vcf \
         -output_dir ${workdir}/ 
        
    """

}
